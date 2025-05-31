import numpy as np
import os.path
import sys


# retorna la version de gaussian
def get_version(file_name):
    version = ""
    with open(file_name, 'r') as file:
        for line in file:
            if "Gaussian(R)" in line:
                line = line.split("Gaussian(R)")[1]
                version = line.strip().split(" ")[0]
                break
    return version


# retorna el metodo de calculo
def get_method(file_name):
    method = ""
    with open(file_name, 'r') as file:
        for line in file:
            if "#" in line.lower() and "irc" in line.lower():
                line = line.strip("\n").split(" ")
                for x in line:
                    if "/" in x:
                        method = x
                        break
                break
    return method


def test_termination(file_name):
    completed = False
    with open(file_name, 'r') as file:
        for line in file:
            if "Normal termination of Gaussian" in line:
                completed = True
    return completed


def IRC_direction(file_name):
    direction = "both"
    with open(file_name, 'r') as file:
        for line in file:
            if "#" in line.lower() and "irc" in line.lower():
                if "forward" in line.lower():
                    direction = "forward"
                elif "reverse" in line.lower():
                    direction = "reverse"
                else:
                    direction = "both"
    return direction


# retorna  la carga y la multiplicidad
def get_ch_mult(file_name):
    ch_mult = "0 1"
    with open(file_name, 'r') as file:
        for line in file:
            if "Multiplicity" in line:
                line = line.strip("\n").split(" ")
                line = [x for x in line if x != ""]
                ch_mult = line[2] + " " + line[5]
                break
    return ch_mult


# retorna una lista con los nombres de los atomos
def atom_list(file_name):
    a_list = []
    with open(file_name, 'r') as file:
        for line in file:
            if "Multiplicity" in line:
                for atomo in file:
                    a = atomo.split(" ")[1]
                    if a == "\n":
                        break
                    else:
                        a_list.append(a)
                break
    return a_list


# retorna el metodo de calculo
def get_recorrect_policy(file_name):
    recorrect = True
    with open(file_name, 'r') as file:
        for line in file:
            if "#" in line.lower() and "irc" in line.lower():
                line = line.strip("\n").split(",")
                for x in line:
                    if "recorrect" in x and "never" in x:
                        recorrect = False
                        break
                break
    return recorrect


# Retorna una lista de las matrices de coordenadas cartesianas de cada punto en Angstroms (A)
# Solo para Gaussian 09
# array(str)
def xyz_matrix_09(file_name, reverse):
    xyz_m = []
    converged = []
    al = atom_list(file_name)
    na = len(al)
    with open(file_name, 'r') as file:
        TS = True
        recorrect = False
        for line in file:
            if TS and "Input orientation:" in line:
                xyz_l = []
                for k in range(0, 4):
                    file.readline()
                for k in range(0, na):
                    line = file.readline()[30:].strip("\n")
                    xyz_l.append(line)
                xyz_m.append(xyz_l)
                converged.append(True)
                TS = False
            elif "Input orientation:" in line and not recorrect:
                xyz_l = []
                for k in range(0, 4):
                    file.readline()
                for k in range(0, na):
                    line = file.readline()[30:].strip("\n")
                    xyz_l.append(line)
                xyz_m.append(xyz_l)
            elif "Delta-x Convergence" in line and not recorrect:
                if "NOT Met" in line:
                    converged.append(False)
                else:
                    converged.append(True)
            elif "Recorrect forced on => Unset convergence flag." in line:
                recorrect = True
            elif "Input orientation:" in line and recorrect:
                xyz_l = []
                for k in range(0, 4):
                    file.readline()
                for k in range(0, na):
                    line = file.readline()[30:].strip("\n")
                    xyz_l.append(line)
                xyz_m[-1] = xyz_l
            elif "Delta-x Convergence" in line and recorrect:
                if "NOT Met" in line:
                    converged[-1] = False
                else:
                    converged[-1] = True
                recorrect = False
            elif "calculation of the REVERSE path" in line:
                xyz_m = list(reversed(xyz_m))
                converged = list(reversed(converged))
            elif "Reaction path calculation complete." in line:
                break
    if not reverse:
        xyz_m = list(reversed(xyz_m))
        converged = list(reversed(converged))
    completed = test_termination(file_name)
    if not completed:
        print("WARNING: Uncompleted calculation")
        xyz_m = xyz_m[:-1]
    xyz_m = np.array(xyz_m)
    xyz_m = xyz_m[np.where(converged)]
    return xyz_m


# Retorna una lista de las matrices de coordenadas cartesianas de cada punto en Angstroms (A)
# Solo para Gaussian 16
# array(str)
def xyz_matrix_16(file_name, reverse):
    xyz_m = []
    al = atom_list(file_name)
    na = len(al)
    with open(file_name, 'r') as file:
        TS = True
        for line in file:
            if TS and "Input orientation:" in line:
                xyz_l = []
                for k in range(0, 4):
                    file.readline()
                for k in range(0, na):
                    line = file.readline()[30:].strip("\n")
                    xyz_l.append(line)
                xyz_m.append(xyz_l)
                TS = False
            elif "Input orientation:" in line:
                xyz_l = []
                for k in range(0, 4):
                    file.readline()
                for k in range(0, na):
                    line = file.readline()[30:].strip("\n")
                    xyz_l.append(line)
                xyz_m.append(xyz_l)
            elif "calculation of the REVERSE path" in line:
                xyz_m = list(reversed(xyz_m))
            elif "Reaction path calculation complete." in line:
                break
    completed = test_termination(file_name)
    if not completed:
        print("WARNING: Uncompleted calculation")
        xyz_m = xyz_m[:-1]
    if not reverse:
        xyz_m = list(reversed(xyz_m))
    xyz_m = np.array(xyz_m)
    return xyz_m


def get_RC(file_name, reverse):
    direction = IRC_direction(file_name)
    RC_list = [0.]
    with open(file_name, 'r') as file:
        for line in file:
            if "NET REACTION COORDINATE UP TO THIS POINT" in line:
                RC_list.append(float(line.split('=')[1].strip()))
            elif "calculation of the REVERSE path" in line:
                RC_list = [-x for x in reversed(RC_list)]
            elif "Reaction path calculation complete." in line:
                break
    RC_list = np.array(RC_list)
    if not reverse and direction != "reverse":
        RC_list = RC_list[::-1] * (-1)
    if len(RC_list) == 0:
        print("Can't find RC coordinates")
        sys.exit(1)
    return RC_list


# def get_RC(file_name, reverse):
#     RC_list = []
#     with open(file_name, 'r') as file:
#         for line in file:
#             if "Summary of reaction path following" in line:
#                 file.readline()
#                 file.readline()
#                 while True:
#                     RC_line = file.readline()
#                     if "---" in RC_line:
#                         break
#                     else:
#                         RC_line = RC_line.strip("\n").split(" ")
#                         RC_line = [x for x in RC_line if x != ""]
#                         RC_list.append(RC_line[2])
#                 break
#     RC_list = list(map(float, RC_list))
#     RC_list = np.array(RC_list)
#     if reverse:
#         RC_list = RC_list[::-1] * (-1)
#     if len(RC_list) == 0:
#         print("Can't find RC coordinates")
#         sys.exit(1)
#     return RC_list


def get_IRC(file_name, reverse=False):
    if not os.path.isfile(file_name):
        print("ERROR: Can't find the file %s" % file_name)
        sys.exit(2)
    else:
        # version = get_version(file_name)
        method = get_method(file_name)
        ch_mult = get_ch_mult(file_name)
        at = atom_list(file_name)
        version = get_version(file_name)
        recorrect = get_recorrect_policy(file_name)
        if version == "09" and recorrect:
          xyz = xyz_matrix_09(file_name, reverse)
        elif version == "16" or not recorrect:
           xyz = xyz_matrix_16(file_name, reverse)
        else:
            print("ERROR: Can't find Gaussian version")
            sys.exit(2)
        RC = get_RC(file_name, reverse)
        return method, ch_mult, at, xyz, RC
