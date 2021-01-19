#! /users/nanje/miniconda3/envs/hlatyping/bin/python3.7
#import pandas library
import pandas as pd
#read in the data
data = pd.read_csv('/cbio/dbs/refpanels/1000G/1000GP_Phase3.sample')
#rename the columns so as to remove spaces
data = data.rename(columns={'ID POP GROUP SEX': 'IDPOPGROUPSEX'})
#split the columns  into 4 different  columns
data[['ID','POP','GROUP','SEX']]= data.IDPOPGROUPSEX.str.split(" ",expand=True,)
#remove the concatenated column
data.pop("IDPOPGROUPSEX")
#extract data from the African samples
AFR = data[data["GROUP"] == "AFR"]
#extract the  ids of the African samples
Afrids=AFR["ID"]
#Write the output to a file
AFR.to_csv('Africansamples.tsv', sep="\t", index=False)
Afrids.to_csv('Africanids.txt', sep=" ", index=False, header=False)


