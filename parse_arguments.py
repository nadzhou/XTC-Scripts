import argparse as ap 


def parse_arguments(parser=None): 
    if parser == None: 
        parser = ap.ArgumentParser()

    parser.add_argument("input_tpr", help="Input gromacs .tpr file")
    parser.add_argument("input_xtc", help="Input gromacs .xtc file")

    parser.add_argument("--input_tpr2", default=None, help="Input gromacs .tpr file")
    parser.add_argument("--input_tpr2", default=None, help="Input gromacs .tpr file")

    parser.add_argument("--input_tpr3", default=None, help="Input gromacs .tpr file")
    parser.add_argument("--input_xtc3", default=None, help="Input gromacs .tpr file")

    parser.add_argument("--input_tpr4", default=None, help="Input gromacs .tpr file")
    parser.add_argument("--input_xtc4", default=None, help="Input gromacs .tpr file")

    parser.add_argument("--input_tpr5", default=None, help="Input gromacs .tpr file")
    parser.add_argument("--input_xtc5", default=None, help="Input gromacs .tpr file")

    parser.add_argument("--input_tpr6", default=None, help="Input gromacs .tpr file")
    parser.add_argument("--input_xtc6", default=None, help="Input gromacs .tpr file")



    return parser.parse_args()


