import os

import numpy as np
from Bio import PDB

def get_single_rna_coordinate(chain_id, pdb_file_path):
    residue_id = ['A', 'U', 'C', 'G', 'DA', 'DU', 'DC', 'DG', 'PSU', 'CBV', '5BU', 'UMS', 'CSL', 'CCC', 'GTP', 'GDP',
                  'A23', 'U37', 'IU']
    rna_coordinate = []
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("MY", pdb_file_path)
    model = structure[0]
    for chain in model:
        if chain.id != chain_id:
            continue
        for residue in chain:
            if residue.get_resname().replace(" ", "") in residue_id:
                index = True
                for atom in residue:
                    if atom.get_name().strip() == "C3'":
                        x_atom, y_atom, z_atom = atom.get_coord()
                        rna_coordinate.append([x_atom, y_atom, z_atom])
                        index = False
                        break
                if index:
                    for atom in residue:
                        if atom.get_name().strip()[0] in ["C", "N", "O", "P"]:
                            x_atom, y_atom, z_atom = atom.get_coord()
                            rna_coordinate.append([x_atom, y_atom, z_atom])
                            index = False
                            break
                if index:
                    print(pdb_file_path)
    return rna_coordinate


def gaussian_kernel(distance, r):
    """Compute Gaussian kernel value."""
    return np.exp(-distance ** 2 / r ** 2)


def compute_laplace_operator(coords, r):
    """
    Compute the discrete Laplace operator for nucleotides.

    Parameters:
        coords: numpy array of shape (N, 3), coordinates of nucleotides.
        r: float, scale factor for the Gaussian kernel.

    Returns:
        laplace_operator: numpy array of shape (N, N), the Laplace operator.
    """
    N = coords.shape[0]
    laplace_operator = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if abs(i - j) > 1:  
                distance = np.linalg.norm(coords[i] - coords[j])
                laplace_operator[i, j] = gaussian_kernel(distance, r)
    
    
    return laplace_operator


def compute_ln(coords, r):
    """
    Compute the Local Nucleotide (LN) descriptor.

    Parameters:
        coords: numpy array of shape (N, 3), coordinates of nucleotides.
        r: float, scale factor for the Gaussian kernel.

    Returns:
        ln_values: numpy array of shape (N,), LN values for each nucleotide.
    """
    laplace_operator = compute_laplace_operator(coords, r)
    ln_values = np.zeros(coords.shape[0])

    for i in range(coords.shape[0]):
        weights = laplace_operator[i]
        weighted_sum = np.sum(weights[:, None] * coords, axis=0)
        if np.sum(weights) == 0:
            weighted_center = coords[i]
        else:
            weighted_center = weighted_sum / np.sum(weights)
        ln_values[i] = np.linalg.norm(coords[i] - weighted_center)

    return ln_values


def encode_nucleotides(coords, quantiles):
    """
    Encode each nucleotide into a vector based on LN at various scales.

    Parameters:
        coords: numpy array of shape (N, 3), coordinates of nucleotides.
        quantiles: list of floats, scale factors for the Gaussian kernel.

    Returns:
        encoded_vectors: numpy array of shape (N, len(quantiles)),
                         encoded feature vectors for each nucleotide.
    """
    encoded_vectors = []

    for r in quantiles:
        ln_values = compute_ln(coords, r)
        encoded_vectors.append(ln_values)

    return np.array(encoded_vectors).T


def get_ln(pdb_file_path, pdb_list):
    quantiles = [np.finfo(float).eps, 0.25, 0.5, 0.75, 1.0]
    all_rna_ln = []
    for pdb in pdb_list:
        
        pdb = pdb[1:]
        chain_specific_path = os.path.join(pdb_file_path, f"{pdb}.pdb")
        four_letter_path = os.path.join(pdb_file_path, f"{pdb[:4]}.pdb")
        resolved_path = chain_specific_path if os.path.exists(chain_specific_path) else four_letter_path
        if not os.path.exists(resolved_path):
            raise FileNotFoundError(f"Missing PDB for {pdb}: tried {chain_specific_path} and {four_letter_path}")
        rna_coordinate = get_single_rna_coordinate(pdb[4], resolved_path)
        coords = np.array(rna_coordinate)
        
        encoded_vectors = encode_nucleotides(coords, quantiles)
        encoded_vectors = np.array(encoded_vectors)
        
        all_rna_ln.append(encoded_vectors)
    
    return all_rna_ln
