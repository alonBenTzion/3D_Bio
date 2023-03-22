import Bio.PDB
from Bio.PDB import PDBParser, PDBList, MMCIFIO

def parse_pdb_file(protein_ID: str, filename: str):
    parser = PDBParser()
    return parser.get_structure(protein_ID, filename)


def get_chain(generator, chain_id):
    for chain in generator:
        if chain.id == chain_id:
            return chain


def get_atoms(chain_object):
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
    return atoms_arr, structure


def align(a_arr_1, a_arr_2, str_1):
    #TODO check what apply does
    #TODO handle warnings

    mmcif_io = MMCIFIO()
    mmcif_io.set_structure(str_1)
    mmcif_io.save("str1.cif")
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(a_arr_1, a_arr_2)
    super_imposer.apply(str_1.get_atoms())
    print(super_imposer.rms)
    mmcif_io.set_structure(str_1)
    mmcif_io.save("aligned_str.cif")
    # io = Bio.PDB.PDBIO()
    # io.set_structure(sample_structure)
    # io.save("1UBQ_aligned.pdb")



if __name__ == '__main__':
    #TODO remember t change to argv1, argv2
    p_1_atoms_array, str_1 = read_pdb('6w4h', 'A')
    p_2_atoms_array, str_2 = read_pdb('6w4h', 'A')
    align(p_1_atoms_array, p_2_atoms_array, str_1)