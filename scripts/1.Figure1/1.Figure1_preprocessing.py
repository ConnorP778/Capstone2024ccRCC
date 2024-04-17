import pandas as pd
import numpy as np
base_dir = '/Users/carternorton/repo/GitHub/Capstone2024ccRCC/'
data = pd.read_csv(f"{base_dir}/Input/TCGA_&_Clinical_Data.tsv", sep='\t', header=0, index_col=0)
#Let's isolate bmi 
bmi = data[['BMI Group',"Sex (TCGA)", "Age at Diagnosis (Clinical)"]]

#Rename Age at Diagnosis
bmi.rename(columns={"Age at Diagnosis (Clinical)":"Age", "Sex (TCGA)" : "TC", "BMI Group":"bmi"}, inplace=True)

#How many samples are missing BMI data?
print(bmi["bmi"].isnull().sum(), "samples are missing BMI data")


#Let's drop those samples
bmi = bmi.dropna()

#Let's select for MRNA data from these patients
mrna = pd.read_csv(f'{base_dir}/Input/TCGA-KIRC.htseq_counts.tsv', sep='\t', header=0, index_col=0)

#Remove any columns that don't end in -01
mrna = mrna[mrna.columns[mrna.columns.str.endswith('-01A')]]

#Now remove the -01 from the column names
mrna.columns = mrna.columns.str.replace('-01A', '')

samples = bmi.index.tolist()

samples = [x for x in samples if x in mrna.columns.tolist()]

print(len(samples), "samples have both BMI and mRNA data")

mrna = mrna[samples]

bmi = bmi.loc[samples]

#Let's clean up the mrna table
#Remove duplicated mrna indices (both the hugo and entrez ids are duplicated it seems)
mrna = mrna[~mrna.index.duplicated(keep='first')]

#Next, remove all rows that are entirely NaN
print("There are", mrna.shape[0], "rows and", mrna.shape[1], "columns in the mrna dataframe")
mrna = mrna.dropna(how='all')
print("There are now", mrna.shape[0], "rows and", mrna.shape[1], "columns in the mrna dataframe")

#Let's save bmi data
#Let's replace index with "-" with "."
bmi.index = bmi.index.str.replace('-', '.')
bmi.to_csv(f'{base_dir}scripts/1.Figure1/Figure1_BMI.csv', sep=',', header=True, index=True)

#Let's save this data
#Remove duplicated mrna indices
gene_count = len(mrna.index.tolist())
mrna = mrna[~mrna.index.duplicated(keep='first')]
print("There were", gene_count - len(mrna.index.tolist()), "duplicated gene indices")

#Remove any index that starts with __
mrna = mrna[~mrna.index.str.startswith('__')]

#this data was log2 (x+1) transformed, so let's reverse that
mrna = np.power(2, mrna) - 1
#Let's convert to integers
mrna = mrna.astype(int)

#Let's save this data
mrna.to_csv(f'{base_dir}/scripts/1.Figure1/KIRC_mRNA_BMI.csv', sep=',', header=True, index=True)
