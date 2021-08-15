#!/users/nanje/miniconda3/bin/python

import sys
sys.path.append('/scratch3/users/nanje/MHC-Imputation-Accuracy/cookHLA/templates')
from MakeGeneticMap.MakeGeneticMap import MakeGeneticMap
from MakeGeneticMap.__main__ import CookHLA_MakeGeneticMap
from src.checkInput import FixInput

import argparse,time
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", default="${prefix}", help="")
parser.add_argument("--hg", default= "18" , help="")
parser.add_argument("--ref", default="${ref}", help="")
parser.add_argument("--output", default="${output}", help="")

args = parser.parse_args()

if __name__ == '__main__':
    CookHLA_MakeGeneticMap(args.prefix, args.hg, args.ref, args.output)
