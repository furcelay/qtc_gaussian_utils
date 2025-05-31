#! /usr/bin/env python

# change the first line to the direction of the python executable

from data_analysis.read_IRC import get_IRC
from data_analysis.read_config import read_config
from data_analysis.write_input import write_one_file, write_multi_files, split_fragments
import argparse
import os
import sys
from shutil import copyfile
import numpy as np


# define a parser to get the sys arguments
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input file name", type=str)
    parser.add_argument("-o", "--output", help="Name of the output file(s)", default="default_IRC",
                        type=str)
    parser.add_argument("-c", "--config", help="Configuration file for input arguments", type=str)
    parser.add_argument("--make_config", help="Make a default configuration file 'default.txt'",
                        choices=["simple", "full"])
    parser.add_argument("-sp", help="Model for single point calculation, Ex: 'hf/6-31+g(2d,p)',"
                                    " the default is the same as the IRC calculation", type=str)
    parser.add_argument("-nbo", help="Do an NBO calculation", action='store_true')
    parser.add_argument("-r", "--reverse", help="Reverse the IRC", action='store_true')
    parser.add_argument("-m", "--multiple", help="Create a single file for each structure",
                        action='store_true')
    parser.add_argument("-p", "--proc", help="Number of processors to use", type=int)
    parser.add_argument("-mem", help="Memory in gigabytes to use", type=int)
    args_ = parser.parse_args()
    return args_


# ====================================================================================
# read arguments
args = parse_arguments()
local = os.path.dirname(os.path.realpath(__file__))

if args.make_config is None:
    # define the calculation
    guess = [False, ""]
    fragments = [False, 1]
    check = [False, ""]
    key = ""
    first_opt = []
    last_opt = []
    title_list = []
    atoms_opt = []
    ch_mult_new = None

    # read de IRC output file
    method, ch_mult, atoms, xyz_list, RC_list = get_IRC(args.input, args.reverse)

    # read de config file
    if args.config is not None:
        name, args.reverse, args.multiple, check, guess, fragments, first_opt, ch_mult_new, \
        atoms_opt, last_opt = read_config(args.config)

        # Check info
        if len(atoms_opt) == 0:
            if fragments[0] and fragments[1] > 1:
                print("Error: Atom options 'Fragment=n' must be specified for fragment splitting")
                sys.exit(1)
            else:
                pass
        elif len(atoms) == len(atoms_opt):
            atoms = atoms_opt
        else:
            print("Error: number of '& Atom options' given not correspond with the number of atoms"
                  " in the IRC")
            sys.exit(1)
        ch_mult = ch_mult_new

    # if not config read inputs
    else:
        name = args.output
        if args.proc is not None:
            first_opt.append("%" + "nproc=%i" % args.proc)
        if args.mem is not None:
            first_opt.append("%" + "mem=%iGB" % args.mem)
        if args.sp is not None:
            method = args.sp
        key += "# " + method
        if args.nbo:
            key += " pop=nboread"
            last_opt.append(" $nbo bndidx $end")
        first_opt.append(key)
        check = [False, name]

    # Names
    for i in range(len(RC_list)):
        title_list.append("Structure number %i at reaction coordinate: %.5f" % (i + 1, RC_list[i]))
    direct = os.getcwd() + "/" + args.output.split(".")[0]

    # Define the Guess for SCF
    if guess[0]:
        guess_list = np.full((len(RC_list), 2), True)
        guess_list[0][1] = False
        if guess[1].lower() == "reactant":
            pass
        elif guess[1].lower() == "product":
            xyz_list = np.flip(xyz_list, axis=0)
            title_list = list(reversed(title_list))
        elif guess[1].lower() == "ts":
            TS = np.where(RC_list == 0)[0][0]
            xyz_list = np.concatenate((np.flip(xyz_list[:TS + 1], axis=0),
                                       xyz_list[TS + 1:]), axis=0)
            title_list = list(reversed(title_list[:TS + 1])) + title_list[TS + 1:]
            guess_list[TS + 1][1] = False
        elif guess[1].split("=")[0].lower() == "n":
            N_gs = int(guess[1].split("=")[1]) - 1
            xyz_list = np.concatenate((np.flip(xyz_list[:N_gs + 1], axis=0),
                                       xyz_list[N_gs + 1:]), axis=0)
            title_list = list(reversed(title_list[:N_gs + 1])) + title_list[N_gs + 1:]
            guess_list[N_gs + 1][1] = False
    else:
        guess_list = np.full((len(RC_list), 2), False)

    # ====================================================================================

    # Write the output
    if not fragments[0]:
        if not args.multiple:
            write_one_file(name, guess_list, first_opt, title_list, ch_mult, atoms, xyz_list, last_opt,
                           check)
        elif args.multiple:
            if not os.path.exists(direct):
                os.makedirs(direct)
                os.chdir(direct)
                write_multi_files(name, guess_list, first_opt, title_list, ch_mult, atoms, xyz_list,
                                  last_opt, check)
            else:
                print("Error: The directory '%s' already exists" % direct)
                sys.exit(1)
    else:
        frag_atoms, frag_xyz, total_ch_mult, frag_ch_mult = split_fragments(fragments[1], atoms, xyz_list, ch_mult)
        write_one_file(name + "_supramolec", guess_list, first_opt, title_list, total_ch_mult, atoms, xyz_list, last_opt,
                       [check[0], check[1] + "_supramolec"])
        for frag in range(fragments[1]):
            write_one_file(name + f"_fragment_{frag+1}", guess_list, first_opt, title_list,
                           frag_ch_mult[frag], frag_atoms[frag], frag_xyz[frag],
                           last_opt, [check[0], check[1] + f"_fragment_{frag+1}"])


# ====================================================================================
elif args.make_config == "simple":
    copyfile(local + "/data_analysis/config-simple.txt", "config-simple.txt")

elif args.make_config == "full":
    copyfile(local + "/data_analysis/config-full.txt", "config-full.txt")
