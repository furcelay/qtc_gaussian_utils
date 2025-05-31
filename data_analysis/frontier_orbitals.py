import re


def get_orbitals_table(file_list):
    orbitals_table = []
    read_table = False
    for line in file_list:
        if 'The electronic state is' in line:
            read_table = True
        elif 'Condensed to atoms' in line:
            break
        elif read_table:
            orbitals_table.append(line)
    return orbitals_table


def split_alpha_beta(orbitals_table):
    alpha = []
    beta = []
    for line in orbitals_table:
        if 'Alpha' in line:
            alpha.append(line)
        elif 'Beta' in line:
            beta.append(line)
    return alpha, beta


def get_occupied_orbitals(file_list):
    occupied_orbitals = []
    for line in file_list:
        orbital_type, energies_str = line.strip().split('--')

        # lista por comprensión
        # re.findall('[0-9]*\\.[0-9]{5}', energies_str):
        # regular expression necesaria si los números no tienen espacios, ej: '558.545632227.27722'
        # lo que hace es buscar un número (caracteres del 0 al 9 -> [0-9]) de cualquier número de dígitos ([0-9]*)
        # seguido de un punto (\\.) y luego un número de 5 dígitos ([0-9]{5})
        # genera una lista con cada aparición de este patrón
        energies = [float(energy) for energy in re.findall('-?[0-9]+\\.[0-9]{5}', energies_str)]

        if 'occ.' in orbital_type:
            occupied_orbitals += energies
    return occupied_orbitals


def get_virtual_orbitals(file_list):
    virtual_orbitals = []
    for line in file_list:
        orbital_type, energies_str = line.strip().split('--')
        energies = [float(energy) for energy in re.findall('-?[0-9]+\\.[0-9]{5}', energies_str)]
        if 'virt.' in orbital_type:
            virtual_orbitals += energies
    return virtual_orbitals


def get_frontier_orbitals(file_list):
    occupied_orbitals = get_occupied_orbitals(file_list)
    virtual_orbitals = get_virtual_orbitals(file_list)
    homo = occupied_orbitals[-1]
    lumo = virtual_orbitals[0]
    return homo, lumo


