import argparse as ap 


def parse_arguments(parser=None): 
    if parser == None: 
        parser = ap.ArgumentParser()

    parser.add_argument("input_tpr", help="Input gromacs .tpr file")
    parser.add_argument("input_xtc", help="Input gromacs .xtc file")

    return parser.parse_args()


