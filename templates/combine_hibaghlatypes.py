#!/usr/bin/env python3

import argparse,time
import csv
# import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--hlatype_files", default="${hlatype_files}", help="Tab separated list of hlatypes")
parser.add_argument("--hlatypes_out", default="${hlatypes_out}", help="")

args = parser.parse_args()

def combine_hibaghlatypes(hlatype_files, hlatypes_out):
    """
    """
    hlatype_files = hlatype_files.split(',')
    outfile = open(hlatypes_out, 'w')
    writer = csv.DictWriter(outfile, fieldnames=["\\tsample.id\\tA.1\\tA.2\\tB.1\\tB.2\\tC.1\\tC.2\\tDQA1.1\\tDQA1.2\\tDQB1.1\\tDQB1.2\\tDRB1.1\\tDRB1.2"])
    writer.writeheader()
    zeros = ['<NA>'] * 6
    fg = 1
    for hla_file in hlatype_files:
        nline = 1
        data = []
        sample = hla_file.split('/')[-1].split('_')[0]
        for line in open(hla_file):
            if nline > 1:
                line = line.strip().split('\\t')
                hla = []
                for it in line[1:7]:
                    it = it.replace('*','')[1:]
                    if it == '':
                        it = '<NA>'
                    hla.append(it)
                data = [str(fg)] + [sample] + hla + zeros
                fg += 1
                outfile.writelines('\\t'.join(data)+'\\n')
            nline += 1
    outfile.close()
    # df = pd.read_csv(hlatypes_out, sep ='\\t')
    #df.insert(0, '', range(1, 1 + len(df)))

    # df = df.reset_index()
    # df.columns[0] = df.index + 1


    # df.to_csv(hlatypes_out, sep = '\\t', index = True)

    # df= pd.read_csv(outfile)
    # df.to_csv(outfile, sep = '\\t', index = True)    

if __name__ == '__main__':
    combine_hibaghlatypes(args.hlatype_files, args.hlatypes_out)