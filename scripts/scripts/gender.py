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

#extractgender of the samples
gender = data["SEX"]
gn = gender.replace({'female': 2})
gy = gn.replace({'male': 1})

gy.to_csv('gender.txt',index=False, header=False)


