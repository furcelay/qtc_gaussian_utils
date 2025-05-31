import re


# function that prints one structure input in a file
def write_single_job(file, first_options, title, ch_mult, atoms, xyz,
                     last_options, guess_i, file_name, check=False, name=""):
    if guess_i[0]:
        print("%chk=" + file_name + ".chk", file=file)
    elif check:
        print("%chk=" + name + ".chk", file=file)
    for line in first_options:
        print(line, file=file)
    if guess_i[1]:
        print("# guess=read", file=file)
    print(" \n %s \n " % title, file=file)
    print(ch_mult, file=file)
    for i in range(len(xyz)):
        print(" %s %s" % (atoms[i], xyz[i]), file=file)
    print(" ", file=file)
    for line in last_options:
        print(line, file=file)
    if len(last_options) != 0:
        print(" ", file=file)


# function that writes one file with all the structures inputs
def write_one_file(file_name, guess_l, first_options, title_list, ch_mult, atoms, xyz_list,
                   last_options, check):
    file = open(file_name + ".gjf", 'w')
    for i in range(len(xyz_list) - 1):
        fname = check[1] + "_" + str(i + 1)
        write_single_job(file, first_options, title_list[i], ch_mult, atoms, xyz_list[i],
                         last_options, guess_l[i], file_name, check[0], fname)
        print("--Link1--", file=file)
    fname = check[1] + "_" + str(len(xyz_list))
    write_single_job(file, first_options, title_list[-1], ch_mult, atoms, xyz_list[-1],
                     last_options, guess_l[-1], file_name, check[0], fname)
    file.close()


# function that writes one input file for each structure
def write_multi_files(file_name, guess_l, first_options, title_list, ch_mult, atoms, xyz_list,
                      last_options, check):
    file_list = open(file_name + "_list.txt", 'w')
    for i in range(len(xyz_list)):
        fname = file_name + "_" + str(i + 1)
        file = open(fname + ".gjf", 'w')
        write_single_job(file, first_options, title_list[i], ch_mult, atoms,
                         xyz_list[i], last_options, guess_l[i], file_name, check[0],
                         check[1] + "_" + str(i + 1))
        print(fname + ".gjf", file=file_list)
        file.close()
    file_list.close()


def split_fragments(n_fragments, atoms_opt, xyz_list, ch_mult):
    frag_index = []
    frag_atoms = [list() for _ in range(n_fragments)]
    frag_xyz_list = [list() for _ in range(n_fragments)]
    for i in range(len(atoms_opt)):
        frag = int(re.findall("Fragment=([1-9]*)", atoms_opt[i])[0]) - 1
        frag_index.append(frag)
        frag_atoms[frag].append(atoms_opt[i])
    for frag in range(n_fragments):
        for j in range(len(xyz_list)):
            frag_xyz = []
            for atom, xyz in enumerate(xyz_list[j]):
                if frag_index[atom] == frag:
                    frag_xyz.append(xyz)
            frag_xyz_list[frag].append(frag_xyz)
    ch_mult = ' '.join(ch_mult.split(',')).split()
    total_ch_mult = ','.join(ch_mult[0:2])
    ch_mult = ch_mult[2:]
    frag_ch_mult = []
    if ch_mult:
        for i in range(0, n_fragments):
            frag_ch_mult.append(','.join(ch_mult[i*2:i*2+2]))
    return frag_atoms, frag_xyz_list, total_ch_mult, frag_ch_mult
