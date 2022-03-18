The folder contains data from the Pathogen decetion database. 
The file "Gene - antibiotic _ Pathogen detection database.xlsx" contains the raw data in first sheet and then step by step cleaned in different sheets. 
From the raw data, only AMR data was selected, the virulance and other stress data like Heat, bile, etc were omitted. The omitted data is still saved in the sheet 7 of the excel file "Gene - antibiotic _ Pathogen detection database.xlsx". 
The allele number like blaCYM_23 were ommited and only the blaCYM remains.Similarly, blaTEM1/2/3 was replaced to blaTEM.  The duplicates were removed. 
If the gene is mentioned to be conferring resistance to more than one antibiotics in one row, the were converted to separate rows. Eg. aph(2'')-I were mentioned to be conferring resitance to Amikacin, Gentamicin and kanamycin. Three different triples were made from it. aph(2'')-I  conferring resitance to antibtiotic (CRA) Amikacin, aph(2'')-I CRA Gentamicin, aph(2'')-I CRA Kanamycin. 
Parenthesis were removed from the gene names, like tet(C) is replaced by tetC. mcr-1.1 were also replaced with mcr-1. 
