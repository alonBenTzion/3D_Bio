import Bio.PDB
from Bio.PDB import PDBParser, PDBList


def parse_pdb_file(protein_ID: str, filename: str):
    parser = PDBParser()
    return parser.get_structure(protein_ID, filename)


def get_chain(generator, chain_id):
    for chain in generator:
        if chain.id == chain_id:
            return chain


def get_atoms(chain_object):
    # TODO check if we need to check 'CA'
    atoms_array = []
    for res in chain_object:
        if 'CA' in res:
            atoms_array.append(res['CA'])
    return atoms_array


def read_pdb(protein_ID: str, chain: str):
    pdb_list = PDBList()
    filename = pdb_list.retrieve_pdb_file(protein_ID, file_format="pdb")
    structure = parse_pdb_file(protein_ID, filename)
    chain_object = get_chain(structure[0].get_chains(), chain)
    atoms_arr = get_atoms(chain_object)
    return atoms_arr


def align(a_arr_1, a_arr_2):
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(a_arr_1, a_arr_2)
    return super_imposer.rms


if __name__ == '__main__':
    p_1_atoms_array = read_pdb('6w4h', 'A')
    p_2_atoms_array = read_pdb('6w4h', 'A')
