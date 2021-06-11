#!/users/nanje/miniconda3/bin/python

import sys
sys.path.append('/scratch3/users/nanje/MHC-Imputation-Accuracy/templates')
from measureAcc.measureAccuracy import CookHLA_measureAcc
import argparse,time
parser = argparse.ArgumentParser()

parser.add_argument("--answer_file", default="${answer_file}", help="")
parser.add_argument("--bglphased", default="${bglphased}", help="")
parser.add_argument("--output", default="${output}", help="")

args = parser.parse_args()

if __name__ == '__main__':
    CookHLA_measureAcc(args.answer_file, args.bglphased, args.output)
