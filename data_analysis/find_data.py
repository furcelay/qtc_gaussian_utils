import numpy as np
from data_analysis.frontier_orbitals import (get_frontier_orbitals, get_orbitals_table, split_alpha_beta,
                                             get_occupied_orbitals, get_virtual_orbitals)


def list_file(file_name):
    try:
        f = open(file_name)
        try:
            f_list = f.readlines()
        finally:
            f.close()
        return f_list
    except (OSError, IOError):
        print("Cant read the file", file_name, "(or file doesn't exists)")


# retorna el metodo de calculo
def method(file_list):
    method = ""
    for line in file_list:
        if "#" in line:
            line = line.strip("\n").split(" ")
            for x in line:
                if "/" in x:
                    method = x.split("/")[0]
                    break
    if method[0:2].upper() == "RO":
        method = method[2:]
    elif (method[0].upper() == "R") or (method[0].upper() == "U"):
        method = method[1:]
    return method


def atom_symbols(file_list):
    a_sym = []
    for k in range(len(file_list)):
        if "Multiplicity" in file_list[k] and "Multiplicity" not in file_list[k+1]:
            # las lineas que siguen contienen la estructura
            for h in range(k+1, len(file_list)):
                a = file_list[h].strip().split(" ")
                if len(a) <= 1:
                    break
                else:
                    a_sym.append(a[0])
            return a_sym


def atom_list(file_list):
    a_list = []
    for i, a in enumerate(atom_symbols(file_list)):
        if a.isalpha():
            # si son simbolos: C -> 1C
            a_list.append(f"{i+1}{a}")
        else:
            # si son numeros atomicos: 6 -> 1(6)
            a_list.append(f"{i+1}({a})")
    return a_list


def energy(file_list, method):
    for line in file_list:
        if "SCF Done" in line:
            e_list = line.split(" ")
            e_list = [x for x in e_list if x != ""]
            if method.upper() in e_list[2].upper():
                return float(e_list[4])    
            else:
                break
    final_block = ""
    block = False
    for line in list(reversed(file_list)):
        if "@" in line:
            block = True
        elif "1\\1\\" in line:
            break
        if block:
            final_block = line.strip("\n").strip(" ") + final_block            
    final_block = final_block.split("\\")
    for x in final_block:
        if (method.upper() + "=") in x:
            E = x.split("=")[1]
            return float(E)
    method = method(file_list).replace("(", "").replace(")", "")
    for x in final_block:
        if (method.upper() + "=") in x:
            E = x.split("=")[1]
            return float(E)


def frontier_orbitals(file_list):
    orbitals_file_list = get_orbitals_table(file_list)

    alpha_list, beta_list = split_alpha_beta(orbitals_file_list)

    homo_alpha, lumo_alpha = get_frontier_orbitals(alpha_list)
    if beta_list:
        homo_beta, lumo_beta = get_frontier_orbitals(beta_list)
    else:
        homo_beta, lumo_beta = homo_alpha, lumo_alpha

    return homo_alpha, lumo_alpha, homo_beta, lumo_beta


def all_orbitals(file_list):
    orbitals_file_list = get_orbitals_table(file_list)

    alpha_list, beta_list = split_alpha_beta(orbitals_file_list)

    alpha_occupied = get_occupied_orbitals(alpha_list)
    alpha_virtual = get_virtual_orbitals(alpha_list)
    if beta_list:
        beta_occupied = get_occupied_orbitals(beta_list)
        beta_virtual = get_virtual_orbitals(beta_list)
    else:
        beta_occupied, beta_virtual = alpha_occupied, alpha_virtual

    return alpha_occupied, alpha_virtual, beta_occupied, beta_virtual


def alpha_orbitals(file_list, n_virtual=None):
    alpha_occupied, alpha_virtual, beta_occupied, beta_virtual = all_orbitals(file_list)
    # concatenated lists
    return alpha_occupied + alpha_virtual[:n_virtual]


def beta_orbitals(file_list, n_virtual=None):
    alpha_occupied, alpha_virtual, beta_occupied, beta_virtual = all_orbitals(file_list)
    # concatenated lists
    return beta_occupied + beta_virtual[:n_virtual]


def mulliken(file_list):
    for k in range(len(file_list)-2):
        if "Mulliken atomic charges" in file_list[k]:
            mull_list = []
            for line in file_list[k+2:]:
                if "Sum of Mulliken" in line:
                    break
                else:
                    at_c = line.strip("\n").split(" ")
                    at_c = [x for x in at_c if x != ""]
                    mull_list.append(at_c[2])
            mull_list = list(map(float, mull_list))
            return mull_list


def nbo_ch(file_list):
    for k in range(len(file_list)-2):
        if "Summary of Natural Population Analysis:" in file_list[k]:
            nbo_list = []
            for line in file_list[k+6:]:
                if "===" in line:
                    break
                else:
                    at_c = line.strip("\n").split(" ")
                    at_c = [x for x in at_c if x != ""]
                    nbo_list.append(at_c[2])
            nbo_list = list(map(float, nbo_list))
            return nbo_list


def Wiberg_BO(file_list, na):
    BO_matrix = []
    for i in range(na):
        BO_matrix.append([])
    N = int(np.ceil(na/9))
    for k in range(len(file_list)):
        if "Wiberg bond index matrix" in file_list[k]:
            n = 0
            i = k
            while n < N:
                n += 1
                i += 4
                for j in range(i, i+na):
                    line = file_list[j]
                    line = line.strip("\n").split(" ")
                    line = [x for x in line if x != ""]
                    BO_matrix[j-i] += line[2:]
                i += na - 1
            BO_list = []
            for i in range(na):
                for j in range(i+1, na):
                    BO_list.append(BO_matrix[i][j])
            BO_list = list(map(float, BO_list))
            return BO_list


def RC(file_list):
    RC = None
    for line in file_list:
        if "reaction coordinate:" in line:
            RC = line.strip("\n").split(": ")[1]
            break
    return float(RC)


def force(file_list, na):
    F_xyz_m = []
    for i, line in enumerate(file_list):
        if "Forces (Hartrees/Bohr)" in line:
            for j in range(i+3, i+3+na):
                line = file_list[j][16:].strip("\n").split(" ")
                line = [x for x in line if x != ""]
                F_xyz_m += line
            return list(map(float, F_xyz_m))


def hessian(file_list, na):
    H_blocks = []
    h_blok = []
    H_matrix = []
    width = 3*na//5 * [5] + [3*na % 5]
    w = 0
    for i in range(len(file_list)):
        if "Force constants in Cartesian coordinates:" in file_list[i]:
            for line in file_list[i+2:]:
                if ("Leave Link  716" in line) or ("Force constants in internal coordinates" in line):
                    while len(h_blok) < 3*na:
                        h_blok.insert(0, width[w]*[0.0])
                    H_blocks.append(h_blok)
                    break
                elif line[6] == " ":
                    while len(h_blok) < 3*na:
                        h_blok.insert(0, width[w]*[0.0])
                    H_blocks.append(h_blok)
                    h_blok = []
                    w += 1
                else:
                    line = line.strip().split(" ")
                    line = [a for a in line if a != ""]
                    line = list(map(lambda x: float(x.replace("D", "E")), line[1:]))
                    while len(line) < width[w]:
                        line.append(0.0)
                    h_blok.append(line)
            break
    if len(H_blocks) > 0:
        for j in range(len(H_blocks[0])):
            h_line = []
            for i in range(len(H_blocks)):
                h_line += H_blocks[i][j]
            H_matrix.append(h_line)
        H_matrix = np.array(H_matrix)
        H_matrix = H_matrix + np.tril(H_matrix, k=-1).T
        return H_matrix.reshape(H_matrix.shape[0]**2)


def vibrations(file_list, na):
    for i in range(len(file_list)):
        if "Harmonic frequencies (cm**-1)" in file_list[i]:
            V = []
            freq = []
            for j, line in enumerate(file_list[i+5:]):
                line = line.strip()
                if line == "":
                    return np.array(freq), np.array(V)
                elif "Frequencies" in line:
                    line = line.split("--")[1].split(" ")
                    line = [float(x) for x in line if x != ""]
                    freq += line
                elif "Atom  AN" in line:
                    vib_bloq = []
                    for vib in file_list[i+5+j+1:i+5+j+1+na]:
                        vib = vib.strip().split(" ")
                        vib = [float(x) for x in vib if x != ""][2:]
                        vib_bloq.append(vib)
                    vib_bloq = np.array(vib_bloq)
                    for k in range(vib_bloq.shape[1]//3):
                        V.append(vib_bloq[:, 3*k:3*(k+1)].reshape(3*na))
            return freq, V


def mass(file_list, na):
    for i in range(len(file_list)):
        if "Thermochemistry" in file_list[i]:
            mass = []
            for line in file_list[i + 3:i + 3 + na]:
                line = line.split("mass")[1].strip()
                mass.append(float(line))
            return np.array(mass)


def coordinates(file_list, na):
    R = []
    for i in range(len(file_list)):
        if "Multiplicity" in file_list[i]:
            for line in file_list[i + 1:i + 1 + na]:
                line = line.strip().split(" ")
                line = [float(x) for x in line[1:] if x != ""]
                R.append(line)
            return np.array(R)


def dipole_moment(file_list):
    for i in range(len(file_list)):
        if "Dipole moment" in file_list[i]:
            md = file_list[i+1].split("Tot=")[1].strip()
            md = float(md)
            return md


def Mayer_BO(file_list, na):
    BO_matrix = []
    for i in range(na):
        BO_matrix.append([])
    N = int(np.ceil(na/6))
    for k in range(len(file_list)):
        if "Mayer Atomic Bond Orders" in file_list[k]:
            n = 0
            i = k
            while n < N:
                n += 1
                i += 2
                for j in range(i, i+na):
                    line = file_list[j]
                    line = line.strip("\n").split(" ")
                    line = [x for x in line if x != ""]
                    BO_matrix[j-i] += line[2:]
                i += na - 1
            BO_list = []
            for i in range(na):
                for j in range(i+1, na):
                    BO_list.append(BO_matrix[i][j])
            BO_list = list(map(float, BO_list))
            return BO_list


def n_fragments(file_list):
    fragments = 1
    for line in file_list:
        if "#" in line:
            line = line.strip("\n").split(" ")
            for x in line:
                if "counterpoise" in x:
                    fragments = int(x.split("=")[1])
            break
    return fragments


def counterpoise(file_list):
    for i in range(len(file_list)):
        if "Counterpoise corrected energy" in file_list[i]:
            counterpoise_corr_line = file_list[i].split(" ")
            counterpoise_corr_line = [x for x in counterpoise_corr_line if x != ""]
            counterpoise_corr = float(counterpoise_corr_line[4])

            bsse_line = file_list[i + 1].split(" ")
            bsse_line = [x for x in bsse_line if x != ""]
            bsse = float(bsse_line[3])

            sum_frag_line = file_list[i + 2].split(" ")
            sum_frag_line = [x for x in sum_frag_line if x != ""]
            sum_frag = float(sum_frag_line[4])

            complx_en_raw_line = file_list[i + 3].split(" ")
            complx_en_raw_line = [x for x in complx_en_raw_line if x != ""]
            # dividimos en 1 kcal/mol para tener todo en hatrees
            complx_en_raw = float(complx_en_raw_line[3]) / 627.5095

            complx_en_corr_line = file_list[i + 4].split(" ")
            complx_en_corr_line = [x for x in complx_en_corr_line if x != ""]
            # dividimos en 1 kcal/mol para tener todo en hatrees
            complx_en_corr = float(complx_en_corr_line[3]) / 627.5095

            return counterpoise_corr, bsse, sum_frag, complx_en_raw, complx_en_corr
