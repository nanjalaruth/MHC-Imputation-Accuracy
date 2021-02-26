#!/usr/bin/env python3

import argparse,time

parser = argparse.ArgumentParser()
parser.add_argument("--hlatype_files", default="${hlatype_files}", help="Comma separated list of hlatypes")
parser.add_argument("--hlatypes_out", default="${hlatypes_out}", help="")

args = parser.parse_args()

def combine_hlatype(hlatype_files, hlatypes_out):
    """
    """
    hlatype_files = hlatype_files.split(',')
    datas = []
    for hla_file in hlatype_files:
        nline = 1
        sample = hla_file.split('_')[0]
        for line in open(hla_file):
            if nline > 1:
                line = line.strip().split('\\t')
                pheno = line[0]
                hla = [ it.replace('*','').replace(':','')[1:] for it in line[1:7] ]
            nline += 1
            data = [ sample, pheno ] + hla
        datas.append('\\t'.join(data))
    print(datas)


if __name__ == '__main__':
    combine_hlatype(args.hlatype_files, args.hlatypes_out)