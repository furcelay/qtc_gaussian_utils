import pandas as pd
import numpy as np
from itertools import combinations, product
from data_analysis import find_data
import os
import sys


def head_combinations(a, sep="-"):
    """
    combinations of an iterable with strings, useful to the header for a property dependent of two atoms like the bond
    order
    example: head_combinations(['C', 'H', 'O']): -> ['C-H', 'C-O', 'O-H']
    :param a:
    :param sep:
    :return:
    """
    comb_list = []
    for comb in list(combinations(a, 2)):
        comb_list.append(sep.join(comb))
    return comb_list


def head_product(a, b, sep="-"):
    """
    product between two iterables with strings, useful to generate the header for a property dependent of
    the coordinates like the force vector
    example: head_product(['C', 'H'], 'xyz'): -> ['C-x', 'C-y', 'C-z', 'H-x', 'H-y', 'H-z']
    :param a:
    :param b:
    :param sep:
    :return:
    """
    comb_list = []
    for comb in list(product(a, b)):
        comb_list.append(sep.join(comb))
    return comb_list


def split_fragments(job_list):
    """
    splits the data of one step in one list for each fragment
    :param job_list: List[str]
    :return: list of lists with the lines for each fragment
    :rtype: List[List[str]]
    """
    fragment_list = [job_list[0]]
    for line in job_list[1:]:
        if "Counterpoise: doing" in line:
            break
        else:
            fragment_list.append(line)
    return fragment_list


class Property:
    """
    This class is a data structure to store a property (energy, coordinates, charges, etc) at each reaction step,
    it uses a function defined at the module find_data and writes the output to a csv file

    Attributes
    ----------
    data_function : Function(List[str]): -> array-like
        function that takes a list of strings corresponding to the lines of the output file of one step in the reaction
        (not the full reaction) and returns an array-like object (float, tuple of floats, array, etc) corresponding to
        this property
    rc_list : array
        a numpy array with the reaction coordinate at each step
    headers : List[str]
        list with the names of the headers for each element of the output, this implicitly defines the dimension of
        this property, if left in None the headers will be auto set numbering the dimensions ['0', '1', ..., 'N']
    is_energy : bool
        True if the property has energy units, is used to scale it to the selected energy unit
    search : bool
        defines if the property should be searched in the file, the default is True but automatically changes to False
        if the property is not found in any step
    extra_args : tuple
        args to be passed to data_function as *args
    data : array
        array containing the value of property at each step

    Methods
    -------
    read(job_list: List[str], rc_idx): -> None
        wrapper for data_function, extracts the property from job_list and assigns it to data[rc_idx]
    write(filename): -> None
        writes the data to a csv file with name filename, uses a pandas data frame for simplicity
    """

    energy_scale = 1

    def __init__(self, data_function, rc_list, headers=None, is_energy=False, search=True, extra_args=None):
        self.data_function = data_function
        self.headers = headers
        self.is_energy = is_energy
        self.search = search
        self.rc_list = rc_list

        if self.headers is None:
            self.data = None
        elif len(self.headers) > 1:
            self.data = np.zeros((len(self.rc_list), len(self.headers)))
        else:
            self.data = np.zeros(len(self.rc_list))

        if extra_args is None:
            self.extra_args = {}
        else:
            self.extra_args = extra_args

    def read(self, job_list, rc_idx):
        value = self.data_function(job_list, **self.extra_args)
        if value is None:
            # stop search for this property
            self.search = False
        elif self.headers is None:
            # auto set the headers from first value, assumes that value has __len__ method
            self.headers = list(range(len(value)))
            self.data = np.zeros((len(self.rc_list), len(self.headers)))

        if self.search:
            self.data[rc_idx] = value
            if self.is_energy:
                self.data[rc_idx] *= self.energy_scale

    def write(self, filename):
        df = pd.DataFrame(data=self.data, index=self.rc_list, columns=self.headers)
        df.index.name = 'RC'
        df.sort_index(inplace=True)
        df.to_csv(filename)


class ReadData:
    """
        This class is a manager composed of Property objects, it reads the full output file and splits it for each step,
        also manages the fragments in counterpoise calculation. To add a new property add an entry in the dictionary
        defined in generate_properties method.

        Attributes
        ----------
        input_file : str
            name of the input file, corresponding to the output of a calculation
        method : str
            energy method to extract (ex: HF) if not given will use the same as the calculation
        reverse : bool
            if True the reaction is inverted
        selected_unit : str
            name of the energy unit to use
        energy_units: Dict{str: float}
             mapping form energy unit name to its conversion from hartrees to se chosen unit
        file_list : List[List[str]]
            list containing a list for each step with the line of the corresponding peace of the file
        rc_list : array
            a numpy array with the reaction coordinate at each step
        atoms : List[str]
            list with the names of the atoms (symbols or numbers depending on the input)
        n_atoms : int
            number of atoms
        n_fragments : int
            number of fragments in a counterpoise calculation

        Methods
        -------
        read_input(): -> None
            read the input file and splits it for each step
        read_rc(): -> None
            reads the reaction coordinate from the input file
        generate_properties(): -> Dict{str: Property}
            defines a dictionary of the properties to be computed, add here the properties you want
        compute_extra_properties(): -> None
            adds new properties to the dictionary that are dependent of other properties
            (ex: derivative of some property)
        read_data_single_fragment(file_list: List[List[str]]): -> Dict{str: Property}
            extracts the properties for each step for a single fragment
        split_job_fragments(): -> Dict{str: Dict{str: Any}}
            splits each job list in file_list in the corresponding fragments in a counterpoise calculation
        read_data_multi_fragments(): -> Dict{str: Dict{str: Property}}
            feeds read_data_single_fragment with the data for each fragment and stores it
        save_data(): -> None
            reads the data and saves the result, its called from __init__
        """

    def __init__(self, input_file, method=None, reverse=False, energy_unit="kcal/mol"):
        self.input_file = input_file
        self.reverse = reverse

        self.selected_unit = energy_unit
        self.energy_units = {"Hartrees": 1, "kcal/mol": 627.5095, "kJ/mol": 2625.4997, "eV": 27.2114}
        Property.energy_scale = self.energy_units[self.selected_unit]

        self.file_list = []
        self.rc_list = np.array([])

        self.read_input()
        self.read_rc()

        self.atoms = find_data.atom_list(self.file_list[0])
        self.n_atoms = len(self.atoms)

        if method is None:
            self.method = find_data.method(self.file_list[0])
        else:
            self.method = method

        self.n_fragments = find_data.n_fragments(self.file_list[0])

        print(f"fragments: {self.n_fragments}, atoms: {self.n_atoms}, method: {self.method}, steps: {len(self.rc_list)}")

        self.save_data()

    def read_input(self):
        job_list = []
        with open(self.input_file, 'r') as file:
            for line in file:
                job_list.append(line)
                if "Normal termination of Gaussian" in line:
                    self.file_list.append(job_list)
                    job_list = []
        if self.reverse:
            self.file_list.reverse()

    def read_rc(self):
        rc_list = []
        for j_list in self.file_list:
            rc_list.append(find_data.RC(j_list))
        self.rc_list = np.array(rc_list)
        if self.reverse:
            self.rc_list *= -1

    def generate_properties(self):
        properties = {
            "energy": Property(find_data.energy, self.rc_list,
                               headers=[f"E"],
                               is_energy=True,
                               extra_args={"method": self.method}),
            "frontier_orbitals": Property(find_data.frontier_orbitals, self.rc_list,
                                          headers=["HOMO-alpha", "LUMO-alpha", "HOMO-beta", "LUMO-beta"],
                                          is_energy=True),
            "alpha_orbitals": Property(find_data.alpha_orbitals, self.rc_list,
                                       headers=None,
                                       is_energy=True,
                                       extra_args={"n_virtual": 10}),
            "beta_orbitals": Property(find_data.beta_orbitals, self.rc_list,
                                      headers=None,
                                      is_energy=True,
                                      extra_args={"n_virtual": 10}),
            "CMull": Property(find_data.mulliken, self.rc_list,
                              headers=self.atoms),
            "CNBO": Property(find_data.nbo_ch, self.rc_list,
                             headers=self.atoms),
            "Wb_BO": Property(find_data.Wiberg_BO, self.rc_list,
                              headers=head_combinations(self.atoms),
                              extra_args={"na": self.n_atoms}),
            "My_BO": Property(find_data.Mayer_BO, self.rc_list,
                              headers=head_combinations(self.atoms),
                              extra_args={"na": self.n_atoms}),
            "force": Property(find_data.force, self.rc_list,
                              headers=head_product(self.atoms, "xyz"),
                              is_energy=True,
                              extra_args={"na": self.n_atoms}),
            "hessian": Property(find_data.hessian, self.rc_list,
                                headers=head_product(head_product(self.atoms, "xyz", sep="_"),
                                                     head_product(self.atoms, "xyz", sep="_")),
                                is_energy=True,
                                extra_args={"na": self.n_atoms}),
            "dipole_moment": Property(find_data.dipole_moment, self.rc_list,
                                      headers=["dm"]),
            "counterpoise": Property(find_data.counterpoise, self.rc_list,
                                     is_energy=True,
                                     headers=["counterpoise_corrected_energy",
                                              "BSSE_energy",
                                              "sum_fragmenst",
                                              "complexation_energy_raw",
                                              "complexation_energy_corrected"])
        }
        return properties

    def compute_extra_properties(self, properties):
        properties["reaction_force"] = Property(None, self.rc_list, headers=["F", "K"],
                                                is_energy=True)
        force = - np.gradient(properties["energy"].data, self.rc_list)
        force_constant = - np.gradient(force, self.rc_list)
        properties["reaction_force"].data[:, 0] = force
        properties["reaction_force"].data[:, 1] = force_constant

        properties["chem_potential"] = Property(None, self.rc_list, headers=["mu", "J"],
                                                is_energy=True)
        mu = properties["frontier_orbitals"].data.sum(axis=1) / 4
        flux = - np.gradient(mu, self.rc_list)
        properties["chem_potential"].data[:, 0] = mu
        properties["chem_potential"].data[:, 1] = flux

        # headers for orbitals
        # homo_alpha, lumo_alpha, homo_beta, lumo_beta = properties["frontier_orbitals"].data[0]
        # alpha_orbitals = properties["alpha_orbitals"].data[0]
        # homo_index = 0
        # for i in range(len(alpha_orbitals)):
        #     if alpha_orbitals[i] == homo_alpha and alpha_orbitals[i + 1] == lumo_alpha:
        #         homo_index = i
        #         break
        # n_occupied = homo_index + 1
        # n_virtual = len(alpha_orbitals) - n_occupied
        # lumo_headers = [f"LUMO+{n}" for n in range(n_virtual)]
        # homo_headers = [f"HOMO-{n}" for n in range(n_occupied - 1, -1, -1)]
        # properties["alpha_orbitals"].headers = homo_headers + lumo_headers

    def read_data_single_fragment(self, file_list):
        properties = self.generate_properties()

        for i, job_list in enumerate(file_list):
            for prop in properties.values():
                if prop.search:
                    prop.read(job_list, i)

        self.compute_extra_properties(properties)
        return properties

    def split_job_fragments(self):
        fragments_jobs = {"total": {"job_list": [],
                                    "key_line": "Counterpoise: doing full-system calculation"}}
        for f in range(self.n_fragments):
            for basis in ["full-system", "fragment"]:
                key_line = f"Counterpoise: doing calculation for fragment   {f + 1} using the basis set of the {basis}"
                fragments_jobs[f"frag_{f + 1}-basis_{basis}"] = {"job_list": [],
                                                                 "key_line": key_line}

        for job_list in self.file_list:
            for i, line in enumerate(job_list):
                for fragment in fragments_jobs.values():
                    if fragment["key_line"] in line:
                        fragment["job_list"].append(split_fragments(job_list[i:]))
        return fragments_jobs

    def read_data_multi_fragments(self):
        fragments_jobs = self.split_job_fragments()
        fragments = {}
        for fragment_name, fragment in fragments_jobs.items():
            fragments[fragment_name] = self.read_data_single_fragment(fragment["job_list"])
        return fragments

    def save_data(self):
        output_dir = os.path.splitext(os.path.split(self.input_file)[1])[0] + "_results"
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            print(f"Directory {output_dir} already exists, data may be overwritten")
            option = input("Continue ([y]/n): ")
            if option.lower() in ["n", "no"]:
                print("Aborted")
                sys.exit(0)

        print("\nReading / Writing ...")

        if self.n_fragments == 1:
            properties = self.read_data_single_fragment(self.file_list)
            for name, prop in properties.items():
                if prop.search:
                    prop.write(os.path.join(output_dir, name + ".csv"))

        else:
            fragments = self.read_data_multi_fragments()
            for frag_name, properties in fragments.items():
                if not os.path.exists(os.path.join(output_dir, frag_name)):
                    os.mkdir(os.path.join(output_dir, frag_name))
                for name, prop in properties.items():
                    if prop.search:
                        prop.write(os.path.join(output_dir, frag_name, name + ".csv"))
            # el archivo de counterpoise queda en el ultimo fragmento, lo movemos al total
            os.rename(os.path.join(output_dir, f"frag_{self.n_fragments}-basis_fragment", "counterpoise.csv"),
                      os.path.join(output_dir, f"total", "counterpoise.csv"))

        with open(os.path.join(output_dir, "info.txt"), 'w', encoding='utf-8') as info_file:
            info_file.write(f"Energy Units: {self.selected_unit} \n")
