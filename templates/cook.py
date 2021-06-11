#!/users/nanje/miniconda3/bin/python

import sys
sys.path.append('/scratch3/users/nanje/MHC-Imputation-Accuracy/templates')
from CookHLA import CookHLA

import argparse,time
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", default="${prefix}", help="")
parser.add_argument("--hg", default= "18" , help="")
parser.add_argument("--ref", default="${ref}", help="")
parser.add_argument("--output", default="${output}", help="")
parser.add_argument("--mem", default="120g", help="")
parser.add_argument("--mp", default="32", help="")

args = parser.parse_args()

if __name__ == '__main__':
    CookHLA(args.prefix, args.hg, args.ref, args.output, args.mem, args.mp)
