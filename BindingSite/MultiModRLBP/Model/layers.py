import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from Rnaglib.config.build_iso_mat import iso_mat


def weighted_MSE(output, target, weight):
    if weight is None:
        return torch.nn.MSELoss()(output, target)
    return torch.mean(weight * (output - target) ** 2)


class RGATLayer(nn.Module):
    def __init__(self, in_feat,
                 out_feat,
                 num_heads=5,
                 num_rels=20,
                 self_loop=True,
                 activation=None,
                 num_bases=-1,
                 iso_mat=iso_mat,
                 return_loss=True,
                 sample_other=0.2,
                 scale_loss=True):
        """

        :param in_feat:
        :param out_feat:
        :param num_rels:
        :param num_heads:
        :param self_loop:
        :param num_bases:
        :param iso_mat:
        :param sample_other: The frequency with which we change the edges in the reconstruction loss
        """
        super(RGATLayer, self).__init__()
        self.num_rels = num_rels
        self.num_heads = num_heads
        self.fc = nn.Linear(in_feat, out_feat, bias=False)
        self.self_interaction = self_loop
        if self_loop:
            self.self_fc = nn.Linear(in_feat, self.num_heads * out_feat, bias=False)
        self.activation = activation
        self.scale_loss = scale_loss
        self.return_loss = return_loss

        
        if num_bases is None or num_bases > self.num_rels:
            self.num_bases = num_rels
            self.use_basis_sharing = False
        else:
            self.num_bases = num_bases
            self.use_basis_sharing = True
        self.attention_weight = nn.Parameter(torch.Tensor(self.num_bases, 2 * out_feat, num_heads))
        if self.use_basis_sharing:
            self.w_comp = nn.Parameter(torch.Tensor(self.num_rels, self.num_bases))

        
        self.sample_other = sample_other
        if sample_other > 0 and iso_mat is None:
            raise ValueError('Cannot use isostericity regularisation if iso_mat is not provided')
        if isinstance(iso_mat, np.ndarray):
            self.iso_mat = torch.from_numpy(iso_mat)

        
        self.reset_parameters()

    def reset_parameters(self):
        """Reinitialize learnable parameters."""
        gain = nn.init.calculate_gain('relu')
        nn.init.xavier_normal_(self.fc.weight, gain=gain)
        if self.self_interaction:
            nn.init.xavier_normal_(self.self_fc.weight, gain=gain)

        
        nn.init.xavier_normal_(self.attention_weight, gain=gain)
        if self.use_basis_sharing:
            nn.init.xavier_normal_(self.w_comp, gain=gain)

    def deactivate_loss(self):
        self.return_loss = False

    def activate_loss(self):
        self.return_loss = True

    def edge_attention(self, edges):
        """
        This is where the relational part hits

        Warning : this needs to stay updated with its sampled copy (could be a factory function
        """
        
        if self.use_basis_sharing:
            attn_fcs = torch.einsum('ab,bcd->acd', self.w_comp, self.attention_weight)
        else:
            attn_fcs = self.attention_weight

        
        z2 = torch.cat([edges.src['z'], edges.dst['z']], dim=1)
        edge_types = edges.data['edge_type']
        w = attn_fcs[edge_types]
        a = torch.bmm(z2.unsqueeze(1), w)
        a = a.squeeze()
        
        return {'attention': F.leaky_relu(a)}

    def sampled_edge_attention(self, edges):
        """
        We need a copy of the function as it causes problem to call a function with arguments
        in the special dgl functions
        """
        
        if self.use_basis_sharing:
            attn_fcs = torch.einsum('ab,bcd->acd', self.w_comp, self.attention_weight)
        else:
            attn_fcs = self.attention_weight

        
        z2 = torch.cat([edges.src['z'], edges.dst['z']], dim=1)

        
        
        
        edge_types = edges.data['edge_type']
        edge_types_perturbation = torch.randint(0, self.num_rels, size=edges.data['edge_type'].shape)
        sampler = torch.distributions.bernoulli.Bernoulli(torch.tensor([self.sample_other]))
        binary_mask = sampler.sample(edge_types.shape).type(torch.bool).squeeze()
        modified_edge_types = edge_types + edge_types_perturbation * binary_mask
        a = np.zeros(10)
        modified_edge_types = torch.remainder(modified_edge_types, torch.tensor(self.num_rels))
        w = attn_fcs[modified_edge_types]
        sampled_attention = torch.bmm(z2.unsqueeze(1), w)
        sampled_attention = sampled_attention.squeeze()
        return {'sampled_attention': F.leaky_relu(sampled_attention), 'modified_edge_types': modified_edge_types}

    def compute_isostericity_loss(self, g):
        """
        Compare the attentions induced distances with the isostericity matrix
        :param g:
        :return:
        """
        native_attentions = g.edata['attention']
        modified_attentions = g.edata['sampled_attention']
        edge_types = g.edata['edge_type']
        modified_edge_types = g.edata['modified_edge_types']

        similarities = self.iso_mat[edge_types, modified_edge_types]

        
        if self.scale_loss:
            with torch.no_grad():
                native_norms = torch.linalg.norm(native_attentions, dim=1)
                modified_norms = torch.linalg.norm(modified_attentions, dim=1)
                scaling = native_norms * modified_norms
                scaling /= scaling.mean()
        else:
            scaling = None
        computed_distances = torch.nn.PairwiseDistance()(native_attentions, modified_attentions)
        computed_similarities = torch.exp(-computed_distances)
        isostericity_loss = weighted_MSE(output=computed_similarities, target=similarities, weight=scaling)
        
        return isostericity_loss

    def message_func(self, edges):
        return {'z': edges.src['z'], 'attention': edges.data['attention']}

    def reduce_func(self, nodes):
        """
        Now use concatenation for the aggregation
        """
        

        all_attentions = nodes.mailbox['attention']  
        all_messages = nodes.mailbox['z']  

        
        scaled_messages = torch.einsum('abi,abj->abij',
                                       all_attentions, all_messages)  
        concatenated = torch.flatten(scaled_messages, start_dim=2)  
        summed_attentions = torch.sum(concatenated, dim=1)  

        if self.self_interaction:
            self_message = nodes.data['self_z']
            
            h = summed_attentions + self_message
        if self.activation:
            h = self.activation(h)
        return {'h': h}

    def forward(self, g, feat):
        
        z = self.fc(feat)
        g.ndata['z'] = z
        if self.self_interaction:
            self_z = self.self_fc(feat)
            g.ndata['self_z'] = self_z

        
        g.apply_edges(self.edge_attention)

        
        if self.return_loss:
            g.apply_edges(self.sampled_edge_attention)
            sample_loss = self.compute_isostericity_loss(g)

        
        g.update_all(self.message_func, self.reduce_func)
        g.ndata.pop('z')
        g.ndata.pop('self_z')
        g.edata.pop('attention')
        if self.return_loss:
            g.edata.pop('sampled_attention')
            g.edata.pop('modified_edge_types')
        h = g.ndata.pop('h')
        if self.return_loss:
            return h, sample_loss
        else:
            return h
