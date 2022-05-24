

import pandas as pd
import sys

# read data
with open('GO_Sal_Raw.txt') as f:
    data = f.readlines()

# process data
predicate = {'C' : 'is part of', 'P' : 'is involved in', 'F': 'has'}
triples = []

for line in data:
    line = line[:-1]
    element = line.split("\t")
    gene = element[2]
    terms = element[3].split('|')
    for term in terms:
        triples.append([gene, predicate[element[4]], term])

# save data
pd_data = pd.DataFrame(triples, columns=['Subject', 'Predicate', 'Object'])
pd_data = pd_data.drop_duplicates()
pd_data = pd_data[pd_data['Object'] != 'molecular_function']
pd_data = pd_data[pd_data['Object'] != 'biological_process']
pd_data = pd_data[pd_data['Object'] != 'cellular_component']

#Exporting data into a csv file
pd_data.to_csv('GO.txt', sep='\t', index=False)
#Exporting data into an Excel file
pd_data.to_excel('GO.xlsx')