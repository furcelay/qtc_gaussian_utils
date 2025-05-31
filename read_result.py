#! /usr/bin/env python

# change the first line to the direction of the python executable

import argparse
from data_analysis.properties import ReadData


# define a parser to get the sys arguments
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input file name", type=str)
    parser.add_argument("-r", "--reverse", help="Reverse the IRC", action='store_true')
    parser.add_argument("-m", "--method",
                        help="Select another energy to extract (Ex: HF)", type=str)
    parser.add_argument("-e", "--energy", help="Energy unit", choices=["Hartrees", "kcal/mol", "kJ/mol", "eV"],
                        default="kcal/mol")
    args_ = parser.parse_args()
    return args_


if __name__ == "__main__":
    args = parse_arguments()
    ReadData(args.input, args.method, args.reverse, args.energy)
    print("Done!")
