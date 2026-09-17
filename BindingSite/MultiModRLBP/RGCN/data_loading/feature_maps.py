"""
Functions to build feature map for each feature type.
"""

import torch


class OneHotEncoder:
    def __init__(self, mapping, num_values=None):
        """
        To one-hot encode this feature.

        :param mapping: This is a dictionnary that gives an index for each possible value.
        :param num_values: If the mapping can be many to one, you should specifiy it here.
        """
        self.mapping = mapping
        self.reverse_mapping = {value: key for key, value in mapping.items()}
        if num_values is None:
            num_values = max(mapping.values())
        self.num_values = num_values

    def encode(self, value):
        """
        Assign encoding of `value` according to known possible values.

        :param value: The value to encode. If missing a default vector of full zeroes is produced.
        """
        x = self.encode_default()
        try:
            ind = self.mapping[value] 
            x[ind] = 1.
            return x
        except KeyError:
            return x

    def encode_default(self):
        x = torch.zeros(self.num_values)
        return x

    def decode(self, one_hot):
        try:
            decoded = self.reverse_mapping[torch.where(one_hot)[0].item()]
            return decoded
        except KeyError:
            return None


class FloatEncoder:

    def __init__(self, default_value=0):
        """
        Utility class to encode floats

        :param default_value: The value to return in case of failure
        """
        self.default_value = default_value

    def encode(self, value):
        """
        Assign encoding of `value` according to known possible values.

        :param value: The value to encode. If missing a default value (by default zero) is produced.
        """
        try:
            return torch.tensor([value], dtype=torch.float)
        except:
            return self.encode_default()

    def encode_default(self):
        return torch.tensor([self.default_value], dtype=torch.float)

    def decode(self, value):
        return value.item()


class BoolEncoder:

    def __init__(self, default_value=False):
        """
        To encode bools. A possible encoding is to have no value in which case it defaults to False.

        :param default_value: To switch the default behavior. This is not recommended because not aligned with the data
        """
        self.default_value = default_value

    def encode(self, value):
        """
        Assign encoding of `value` according to known possible values.

        :param value: The value to encode. If missing the default value (False by default) is produced.
        """
        if value is None:
            return self.encode_default()
        
        if not isinstance(value, bool):
            x = torch.tensor([True], dtype=torch.float)
            return x
        x = torch.tensor([value], dtype=torch.float)
        return x

    def encode_default(self):
        x = torch.tensor([self.default_value], dtype=torch.float)
        return x

    def decode(self, value):
        return value.item()


class ListEncoder:
    def __init__(self, list_length):
        """
        To encode lists, cast them as tensor if possible, otherwise just return zeroes.

        :param list_length: We need the lists to be fixed length
        """
        size = [list_length]
        self.default_value = torch.zeros(size=size, dtype=torch.float)

    def encode(self, value):

        """
        Assign encoding of `value` according to known possible values.

        :param value: The value to encode. If missing the default value (A list of zeros) is produced.
        """
        if value is None or any([val is None for val in value]):
            return self.encode_default()
        else:
            try:
                x = torch.tensor(value, dtype=torch.float)
            except:
                return self.encode_default()
        return x

    def encode_default(self):
        return self.default_value

    def decode(self, value):
        return value.item()















NODE_FEATURE_MAP = {
    'index': None,  
    'index_chain': None,  
    'chain_name': None,  
    'nt_resnum': None,  
    "nt_name": None,  
    'nt_code': OneHotEncoder(mapping={'A': 0, 'U': 1, 'C': 2, 'G': 3, 'a': 0, 'u': 1, 'c': 2, 'g': 3}, num_values=4),
    "nt_id": None,  
    "nt_type": None,  
    "dbn": OneHotEncoder(mapping={'(': 0, '.': 1, ')': 2, '[': 3, ']': 4, '{': 5, '}': 6,
                                  '>': 7, '<': 8, 'A': 9, 'a': 9, 'B': 10, 'b': 10}),
    "summary": None,  
    "alpha": FloatEncoder(),
    "beta": FloatEncoder(),
    "gamma": FloatEncoder(),
    "delta": FloatEncoder(),
    "epsilon": FloatEncoder(),
    "zeta": FloatEncoder(),
    "epsilon_zeta": FloatEncoder(),
    "bb_type": OneHotEncoder(mapping={'BI': 0, '--': 1, 'BII': 2}),
    "chi": FloatEncoder(),
    "glyco_bond": OneHotEncoder(mapping={'anti': 0, '--': 1, 'syn': 2}),
    "C5prime_xyz": ListEncoder(list_length=3),
    "P_xyz": ListEncoder(list_length=3),
    
    "form": OneHotEncoder(mapping={'anti': 0, '--': 1, 'syn': 2}),
    "ssZp": FloatEncoder(),
    "Dp": FloatEncoder(),
    "splay_angle": FloatEncoder(),
    "splay_distance": FloatEncoder(),
    "splay_ratio": FloatEncoder(),
    "eta": FloatEncoder(),
    "theta": FloatEncoder(),
    "eta_prime": FloatEncoder(),
    "theta_prime": FloatEncoder(),
    "eta_base": FloatEncoder(),
    "theta_base": FloatEncoder(),
    "v0": FloatEncoder(),
    "v1": FloatEncoder(),
    "v2": FloatEncoder(),
    "v3": FloatEncoder(),
    "v4": FloatEncoder(),
    "TotalAsa": FloatEncoder(),
    "PolarAsa": FloatEncoder(),
    "ApolarAsa": FloatEncoder(),  
    "amplitude": FloatEncoder(),
    "phase_angle": FloatEncoder(),
    
    "puckering": OneHotEncoder(
        mapping={"C3'-endo": 0, "C4'-exo": 1, "C2'-exo": 2, "C3'-exo": 3, "O4'-endo": 4, "O4'-exo": 5,
                 "C4'-endo": 6, "C2'-endo": 7, "C1'-exo": 8, "C1'-endo": 9, '': 10,'fixed':11}),
    "sugar_class": OneHotEncoder(mapping={"~C3'-endo": 0, "~C2'-endo": 1, '--': 2,'fixed':3}),
    "bin": OneHotEncoder(mapping={'inc': 0, '33p': 1, 'trig': 2, '32p': 3, '22t': 4, '23p': 5, '33t': 6, '32t': 7,
                                  '23m': 8, '23t': 9, '22p': 10, '22m': 11, '33m': 12, '32m': 13,'fixed':14}),
    "cluster": OneHotEncoder(mapping={'__': 0, '1a': 1, '1L': 2, '!!': 3, '1[': 4, '0a': 5, '1c': 6, '&a': 7, '1e': 8,
                                      '1g': 9, '9a': 10, '7a': 11, '1b': 12, '2a': 13, '0b': 14, '4d': 15, '6g': 16,
                                      '4b': 17, '6n': 18, '5n': 19, '1m': 20, '1z': 21, '2[': 22, '3d': 23, '5j': 24,
                                      '6j': 25, '1t': 26, '2g': 27, '7d': 28, '2h': 29, '6d': 30, '7p': 31, '2o': 32,
                                      '2u': 33, '1o': 34, '2z': 35, '5z': 36, '6p': 37, '8d': 38, '3a': 39, '1f': 40,
                                      '#a': 41, '3b': 42, '4n': 43, '5d': 44, '0i': 45, '4a': 46, '7r': 47, '5p': 48,
                                      '4p': 49, '4g': 50, '5q': 51, '5r': 52, '0k': 53, '4s': 54}),
    "suiteness": FloatEncoder(),
    "filter_rmsd": FloatEncoder(),
    "frame_rmsd": FloatEncoder(),
    "frame_origin": ListEncoder(list_length=3),
    "frame_x_axis": ListEncoder(list_length=3),
    "frame_y_axis": ListEncoder(list_length=3),
    "frame_z_axis": ListEncoder(list_length=3),
    "frame_quaternion": ListEncoder(list_length=4),
    "sse_sse": None,  
    "binding_protein": BoolEncoder(),
    "binding_ion": BoolEncoder(),
    "binding_small-molecule": BoolEncoder(),
    
    
    
    
    "binding_protein_id": None,  
    "binding_protein_nt-aa": None,  
    "binding_protein_nt": None,  
    "binding_protein_aa": None,  
    "binding_protein_Tdst": FloatEncoder(),
    "binding_protein_Rdst": FloatEncoder(),
    "binding_protein_Tx": FloatEncoder(),
    "binding_protein_Ty": FloatEncoder(),
    "binding_protein_Tz": FloatEncoder(),
    "binding_protein_Rx": FloatEncoder(),
    "binding_protein_Ry": FloatEncoder(),
    "binding_protein_Rz": FloatEncoder(),
    "is_modified": BoolEncoder(),
    "is_broken": BoolEncoder(),
}


EDGE_FEATURE_MAP = {
    "LW": None,  
    "backbone": BoolEncoder(),
    "nt1": None,  
    "nt2": None,  
    "bp": None,  
    "name": None,  
    "Saenger": None,  
    "DSSR": None,  
}


def build_node_feature_parser(asked_features=None, node_feature_map=NODE_FEATURE_MAP):
    """
    This function will load the predefined feature maps available globally.
    Then for each of the features in 'asked feature', it will return an encoder object for each of the asked features
    in the form of a dict {asked_feature : EncoderObject}

    If some keys don't exist, will raise an Error. However if some keys are present but problematic,
    this will just cause a printing of the problematic keys
    :param asked_features: A list of string keys that are present in the encoder
    :return: A dict {asked_feature : EncoderObject}
    """
    if asked_features is None:
        return {}
    if any([feature not in node_feature_map for feature in asked_features]):
        problematic_keys = tuple([feature for feature in asked_features if feature not in node_feature_map])
        raise ValueError(f'{problematic_keys} were asked as a feature or target but do not exist')

    
    encoding_features = [feature for feature in asked_features if node_feature_map[feature] is not None]
    if len(encoding_features) < len(asked_features):
        unencodable_keys = [feature for feature in asked_features if node_feature_map[feature] is None]
        print(f'{unencodable_keys} were asked as a feature or target but do not exist')
    subset_dict = {k: node_feature_map[k] for k in encoding_features}
    return subset_dict


def build_edge_feature_parser(asked_features=None):
    raise NotImplementedError
