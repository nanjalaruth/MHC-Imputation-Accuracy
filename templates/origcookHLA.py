#!/opt/conda/bin/python3.6

import argparse,time 
parser = argparse.ArgumentParser()
parser.add_argument("--answer_file", default="${answer_file}", help="")
parser.add_argument("--bglphased", default="${bglphased}", help="")
parser.add_argument("--output", default="${output}", help="")

args = parser.parse_args()

def measure_accuracy(answer_file, bglphased, output):
    """
    """
    cd /usr/local/bin
    python -m measureAcc answer_file bglphased output 

if __name__ == '__main__':
    measure_accuracy(args.answer_file, args.bglphased, args.output)
