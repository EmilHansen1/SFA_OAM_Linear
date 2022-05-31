import numpy as np
from extra_packages.OutputInterface import OutputInterface

def Rz(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0.],
                     [np.sin(theta), np.cos(theta), 0.],
                     [0., 0., 1.]])

def Ry(theta):
    return np.array([[np.cos(theta), 0., np.sin(theta)],
                     [0., 1., 0.],
                     [-np.sin(theta), 0., np.cos(theta)]])

def align_along_z(atom1, atom2, atom_pos, center_between=False):
    """
    Align to of the atoms along the z direction. atom2 will be pointing along +z and atom1 in (0,0,0).
    If center_between, then output will be centered between atom1 and atom2 instead.
    """
    atom_pos = atom_pos.copy()  # Don't mess with the original one!
    # Remember indecies should be shifted down one compared to the dictionary from atom_info()
    atom1_i = atom1 - 1
    atom2_i = atom2 - 1

    # First shift atom1 to center
    atom_pos = atom_pos - atom_pos[atom1_i]

    # Find angles and align along z-axis
    atom2_pos = atom_pos[atom2_i]
    theta = np.arccos(atom2_pos[2] / np.linalg.norm(atom2_pos))
    phi = np.arctan2(atom2_pos[1], atom2_pos[0])
    atom_pos = Rz(-phi) @ atom_pos.T
    atom_pos = (Ry(-theta) @ atom_pos).T

    # Recenter if wanted
    if center_between:
        atom_pos = atom_pos - atom_pos[atom2_i]/2

    return atom_pos

def align_along_x(atom1, atom2, atom_pos):
    """
    Rotate atom1 and atom2 along x axis by rotating around z axis. atom2 will be at +x and atom2 at -x.
    """
    atom_pos = atom_pos.copy()  # Don't mess with the original one!
    # Remember indecies should be shifted down one compared to the dictionary from atom_info()
    atom1_i = atom1 - 1
    atom2_i = atom2 - 1

    # Find phi and twist around z
    diff_vec = atom_pos[atom2_i] - atom_pos[atom1_i]
    phi = np.arctan2(diff_vec[1], diff_vec[0])

    return (Rz(-phi) @ atom_pos.T).T

def full_align(atom1_z, atom2_z, atom1_y, atom2_y, atom_pos, center_between=False):
    """
    Align molecule first along z direction and then twist into the xz plane
    """
    atom_pos = atom_pos.copy()
    atom_pos = align_along_z(atom1_z, atom2_z, atom_pos, center_between=center_between)
    atom_pos = align_along_x(atom1_y, atom2_y, atom_pos)
    return atom_pos


def find_closest(target_atom, atom_pos, element_dict, filter_type='None'): 
    """
    Finds the atoms closest to choosen atom.
    """
    target_i = target_atom - 1 
    element_list = np.array([element_dict[key][0] for key in element_dict])
    atom_nr = np.array(list(element_dict.keys()))
    dist_array = np.linalg.norm(atom_pos - atom_pos[target_i], axis=1)
    sort_mask = np.argsort(dist_array)
    print(f'\nAtoms closest to atom nr. {atom_nr[target_i]} ({element_list[target_i]}):')
    print(atom_nr[sort_mask])
    #print(sort_mask)
    #print(element_list[sort_mask])
    
def shift_to_geo_center(target_atoms, atom_pos): 
    """
    Shifts the center of the molecule to the geometric center of the atoms given in the list target_atoms. 
    """
    target_atoms = np.array(target_atoms) - 1
    geo_center = np.sum(atom_pos[target_atoms], axis=0) / 2
    return atom_pos - geo_center
    
def print_input_format(atom_pos, element_dict):
    """
    Print atom information in a form suitable for GAMESS input file
    """
    print('\nNew configuartion ready to put in GAMESS .inp file:\n')
    for i, pos in enumerate(atom_pos):
        print(f'{element_dict[i+1][0]} \t {element_dict[i+1][1]} \t {pos[0]:.10f} \t {pos[1]:.10f} \t {pos[2]:.10f}')


file_name = 'output_files/norcamphor.out'
output = OutputInterface(file_name, convert_to_bohr=False)
output.atom_info()
pos_dict = output.position
atom_pos = []

for i in output.position:
    atom_pos.append(pos_dict[i])
atom_pos = np.array(atom_pos)

find_closest(12, atom_pos, output.elements)
atom_pos = align_along_z(1, 12, atom_pos, center_between=True)
atom_pos = align_along_x(2, 3, atom_pos)

atom_pos = shift_to_geo_center([2,3], atom_pos)

print_input_format(atom_pos, output.elements)
#print(align_xz_plane(3,1,atom_pos))
#print_input_format(atom_pos, output.elements)
