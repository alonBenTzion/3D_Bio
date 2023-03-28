import Bio.PDB as PDB
from Bio.PDB import PDBParser, PDBList, MMCIFIO
import sys


def parse_pdb_file(protein_ID: str, filename: str):
    parser = PDBParser()
    return parser.get_structure(protein_ID, filename)


def get_chain(generator, chain_id):
    """
    Get a specific chain
    """
    for chain in generator:
        if chain.id == chain_id:
            return chain


def get_atoms(chain_object):
    """
    Get array of 'CA' atoms from a given chain
    """
    atoms_array = []
    for res in chain_object:
        if 'CA' in res:
            atoms_array.append(res['CA'])
    return atoms_array


def preprocess_protein(protein_ID: str, chain: str):
    """
    Read, parse and get list of atoms of a given protein
    """
    pdb_list = PDBList()
    filename = pdb_list.retrieve_pdb_file(protein_ID, file_format="pdb")
    structure = parse_pdb_file(protein_ID, filename)
    chain_object = get_chain(structure[0].get_chains(), chain)
    atoms_arr = get_atoms(chain_object)
    return atoms_arr, structure


def align(a_arr_1, a_arr_2, str_1, str_2):
    """
    Align two proteins and save each protein as cif file
    """
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(a_arr_1, a_arr_2)
    super_imposer.apply(a_arr_2)
    io = MMCIFIO()
    io.set_structure(str_1)
    io.save(f"{str_1.id}.cif")
    io.set_structure(str_2)
    io.save(f"{str_2.id}.cif")
    print(super_imposer.rms)


if __name__ == '__main__':
    p_1_atoms_array, str_1 = preprocess_protein(sys.argv[1], sys.argv[2])
    p_2_atoms_array, str_2 = preprocess_protein(sys.argv[3], sys.argv[4])
    align(p_1_atoms_array, p_2_atoms_array, str_1, str_2)
