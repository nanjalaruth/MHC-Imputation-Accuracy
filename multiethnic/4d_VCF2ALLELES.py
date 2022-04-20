#-*- coding: utf-8 -*-
import os, sys, re
import numpy as np
import pandas as pd

def toAlleles(_input, _out_prefix):

    ### Loading input ###

    df_input = pd.read_csv(_input, sep='\s+', header=0, dtype=str)
    # print("df_input:\n{}\n".format(df_input))

    df_HLA = df_input.iloc[:, 0].str.extract(r'^HLA_(\w+)\*(.+)$', expand=True)
    df_HLA.columns = ['HLA', 'Allele']
    # print("df_HLA:\n{}\n".format(df_HLA))

    # pattern for 1-field.
    p_1field = re.compile(r'^\d{2,3}') 

    ### Main iteration ###

    count = 0
    l_RETURN = [] # for each sample
    for SampleID, sr_col in df_input.iloc[:, 1:].iteritems():

        df_chrs = sr_col.str.split('|', expand=True)
        df_chrs.columns = ['chr1', 'chr2']

        df_temp = pd.concat([df_HLA, df_chrs], axis=1)

        count1 = 0
        for hla, df in df_temp.groupby('HLA'):
            # print("HLA: {}".format(hla))
            f_chr1 = df['chr1'].values == '1'
            f_chr2 = df['chr2'].values == '1'

            arr_Allele = df['Allele'].values
           

            allele1 = ';'.join(arr_Allele[f_chr1]) if np.any(f_chr1) else ''
            allele2 = ';'.join(arr_Allele[f_chr2]) if np.any(f_chr2) else ''

            allele1_1field = p_1field.match(allele1).group(0) if bool(p_1field.match(allele1)) else ''
            allele2_1field = p_1field.match(allele2).group(0) if bool(p_1field.match(allele2)) else ''

            allele1_2field = allele1.split(';')[1] if len(allele1.split(';')) > 2 else ' '
            allele1_2field = allele1_2field.replace(':','')
            allele2_2field = allele2.split(';')[1] if len(allele2.split(';')) > 2 else ' '
            allele2_2field  = allele2_2field.replace(':','')
           
           
            l_temp = [0, SampleID, hla, ','.join([allele1_1field, allele2_1field]), ','.join([allele1_2field, allele2_2field])]
            
            l_RETURN.append(l_temp)

            count1 +=1
            # if count1 >= 1: break

        count += 1
        # if count >= 1: break


    df_RETURN = pd.DataFrame(l_RETURN)
    print("df_RETURN:\n{}\n".format(df_RETURN))

    df_RETURN.to_csv(_out_prefix+'.alleles', sep='\t', header=False, index=False)
    print(_out_prefix+'.alleles')

    return _out_prefix+'.alleles'



if __name__ == '__main__':

    [_input, _out] = sys.argv[1:]

    toAlleles(_input, _out)
