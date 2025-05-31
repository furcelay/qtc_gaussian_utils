import sys


def read_config(config_file):
    name = ""
    reverse = False
    multiple = False
    checkpoint = False
    check_name = ""
    read_guess = False
    guess_start = ""
    split_frag = False
    n_frag = 1
    link0 = []
    keywords = []
    charge_mult = ""
    atom_opt = []
    last_opt = []
    with open(config_file, 'r') as file:
        config = file.readlines()
    for i, line in enumerate(config):
        if line[0] == "&":
            if "name" in line:
                name = config[i + 1].strip("\n")
            elif "reverse" in line:
                reverse = config[i + 1].strip().lower() == "true"
            elif "multiple" in line:
                multiple = config[i + 1].strip().lower() == "true"
            elif "checkpoint" in line:
                checkpoint = config[i + 1].strip().lower() == "true"
                if checkpoint:
                    check_name = config[i + 2].strip("\n")
            elif "Read Guess" in line:
                read_guess = config[i + 1].strip().lower() == "true"
                if read_guess:
                    guess_start = config[i + 2].strip("\n")
            elif "Split Fragments" in line:
                split_frag = config[i + 1].strip().lower() == "true"
                if split_frag:
                    n_frag = int(config[i + 2].strip("\n"))
            elif "Link 0" in line:
                for link in config[i+1:]:
                    if (link[0] == "!") or (link[0] == "&") or (link.strip("\n").strip(" ") == ""):
                        break
                    else:
                        link0.append(link.strip("\n"))
            elif "Keywords" in line:
                for key in config[i+1:]:
                    if (key[0] == "!") or (key[0] == "&") or (key.strip("\n").strip(" ") == ""):
                        break
                    else:
                        keywords.append(key.strip("\n"))
            elif "Charge and Multiplicity" in line:
                charge_mult = config[i+1].strip("\n")
            elif "Atom options" in line:
                for at in config[i+1:]:
                    if (at[0] == "!") or (at[0] == "&") or (at.strip("\n").strip(" ") == ""):
                        break
                    else:
                        atom_opt.append(at.strip("\n"))
            elif "After coordinates options" in line:
                for lo in config[i+1:]:
                    if (lo[0] == "!") or (lo[0] == "&") or (lo.strip("\n").strip(" ") == ""):
                        break
                    else:
                        last_opt.append(lo.strip("\n"))
    if name == "":
        print("'$ name' must be specified")
        sys.exit(1)
    elif len(keywords) == 0:
        print("'$ keywords' must be specified")
        sys.exit(1)
    elif charge_mult == "":
        print("'$ Charge and Multiplicity' must be specified")
        sys.exit(1)
    if check_name == "" or check_name.lower() == "default":
        check_name = name
    first_opt = link0 + keywords
    return (name, reverse, multiple, [checkpoint, check_name], [read_guess, guess_start], [split_frag, n_frag],
            first_opt, charge_mult, atom_opt, last_opt)
