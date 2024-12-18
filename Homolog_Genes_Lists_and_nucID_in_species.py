#Get homolog genes and nucleotide IDs for various species for human genes which have AUG dORFs in the 3' UTR

#import required module
import xml.etree.ElementTree as ET
import pandas as pd #Import pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #prevent FutureWarning note printed about df.append for each query

#import .csv dataframe and read into df then take gene name column and convert into list with duplicates removed
df_RefSeq_Three_prime = pd.read_csv("Human_RefSeq_dORF_3UTRs.csv") #convert human accession csv file in current directory to dataframe
HumanGenes = df_RefSeq_Three_prime['Gene Name'].to_list()
HumanGenes = list(set(HumanGenes))

print('Human Gene list created to search Homologene Database')

#Create empty lists ready to store genes
Hsapiens_homolog_genes = list()
Ptroglodytes_homolog_genes = list()
Mmulatta_homolog_genes = list()
Clupus_homolog_genes = list()
Btaurus_homolog_genes = list() 
Mmusculus_homolog_genes = list()
Rnorvegicus_homolog_genes = list()
Ggallus_homolog_genes = list()
Xtropicalis_homolog_genes = list()
Drerio_homolog_genes = list()
Dmelanogaster_homolog_genes = list()
Agambiae_homolog_genes = list()
Celegans_homolog_genes = list()
Scerevisiae_homolog_genes = list()
Klactis_homolog_genes = list()
Egossypii_homolog_genes = list()
Spombe_homolog_genes = list()
Moryzae_homolog_genes = list()
Ncrassa_homolog_genes = list()
Athaliana_homolog_genes = list()
Osativa_homolog_genes = list()
#parse xml file
HomologeneXML = ET.parse('Homologene_Database.xml')
#get first element from xml file
Homologene_Entries = HomologeneXML.getroot()
#loop through each separate entry in Homologene Database
for HGEntry in Homologene_Entries.iter('HG-Entry'):
    #Loop through each set of genes listed under each entry
    for HGGenes in HGEntry.iter('HG-Entry_genes'):
        #Convert xml data for each genes list in each entry into string format
        HGGenes_string = ET.tostring(HGGenes, encoding='utf8').decode('utf8')
        #Create loop to only carry out further function if the entry genes include human gene using taxID 9606
        if '>9606<' in HGGenes_string:
            #loop through each gene in gene list
            for HGGene in HGGenes.iter('HG-Gene'):
                #collect gene symbol for each gene in list
                for HGGeneSymbol in HGGene.iter('HG-Gene_symbol'):
                    gene_symbol = HGGeneSymbol.text
                #collect nucleotide ID for each gene in list
                for HGGeneNucAcc in HGGene.iter('HG-Gene_nuc-acc'):
                    nucleotide_id = HGGeneNucAcc.text
                #collect species taxID for each gene in list
                for HGSpecies in HGGene.iter('HG-Gene_taxid'):
                    species_taxID = HGSpecies.text
                    #if taxID is human add to human genes list
                    if species_taxID == '9606':
                        Hsapiens_homolog_genes.append(gene_symbol)
                    #if taxID is not human add gene to homologs list with human homolog gene by taking the last human gene to be found which will be from this entry
                    elif species_taxID == '9598':
                        Ptroglodytes_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '9544':
                        Mmulatta_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '9615':
                        Clupus_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '9913':
                        Btaurus_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '10090':
                        Mmusculus_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '10116':
                        Rnorvegicus_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '9031':
                        Ggallus_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '8364':
                        Xtropicalis_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '7955':
                        Drerio_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '7227':
                        Dmelanogaster_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '7165':
                        Agambiae_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '6239':
                        Celegans_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '4932':
                        Scerevisiae_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '28985':
                        Klactis_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '33169':
                        Egossypii_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '4896':
                        Spombe_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '318829':
                        Moryzae_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '5141':
                        Ncrassa_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '3702':
                        Athaliana_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    elif species_taxID == '4530':
                        Osativa_homolog_genes.append([Hsapiens_homolog_genes[-1], gene_symbol, nucleotide_id])
                    else:
                        continue
                    #show progress of script
                    print('\r\033[K', end='')
                    print("Homologene Database Entries Searched -", len(Hsapiens_homolog_genes), end='\r')
        else:
            continue

print('Homologene Database Searching completed')
#create lists for filtered homolog genes
F_Hsapiens_homolog_genes = list()
F_Ptroglodytes_homolog_genes = list()
F_Mmulatta_homolog_genes = list()
F_Clupus_homolog_genes = list()
F_Btaurus_homolog_genes = list() 
F_Mmusculus_homolog_genes = list()
F_Rnorvegicus_homolog_genes = list()
F_Ggallus_homolog_genes = list()
F_Xtropicalis_homolog_genes = list()
F_Drerio_homolog_genes = list()
F_Dmelanogaster_homolog_genes = list()
F_Agambiae_homolog_genes = list()
F_Celegans_homolog_genes = list()
F_Scerevisiae_homolog_genes = list()
F_Klactis_homolog_genes = list()
F_Egossypii_homolog_genes = list()
F_Spombe_homolog_genes = list()
F_Moryzae_homolog_genes = list()
F_Ncrassa_homolog_genes = list()
F_Athaliana_homolog_genes = list()
F_Osativa_homolog_genes = list()

#from previously compiles lists only keep them if the human gene homolog matches the list of human genes originally provided
for x in Hsapiens_homolog_genes:
    if x in HumanGenes:
        F_Hsapiens_homolog_genes.append(x)
    else:
        continue

for y in Ptroglodytes_homolog_genes:
    if y[0] in HumanGenes:
        F_Ptroglodytes_homolog_genes.append(y)
    else:
        continue

for y in Mmulatta_homolog_genes:
    if y[0] in HumanGenes:
        F_Mmulatta_homolog_genes.append(y)
    else:
        continue

for y in Clupus_homolog_genes:
    if y[0] in HumanGenes:
        F_Clupus_homolog_genes.append(y)
    else:
        continue

for y in Btaurus_homolog_genes:
    if y[0] in HumanGenes:
        F_Btaurus_homolog_genes.append(y)
    else:
        continue

for y in Mmusculus_homolog_genes:
    if y[0] in HumanGenes:
        F_Mmusculus_homolog_genes.append(y)
    else:
        continue

for y in Rnorvegicus_homolog_genes:
    if y[0] in HumanGenes:
        F_Rnorvegicus_homolog_genes.append(y)
    else:
        continue

for y in Ggallus_homolog_genes:
    if y[0] in HumanGenes:
        F_Ggallus_homolog_genes.append(y)
    else:
        continue

for y in Xtropicalis_homolog_genes:
    if y[0] in HumanGenes:
        F_Xtropicalis_homolog_genes.append(y)
    else:
        continue

for y in Drerio_homolog_genes:
    if y[0] in HumanGenes:
        F_Drerio_homolog_genes.append(y)
    else:
        continue

for y in Dmelanogaster_homolog_genes:
    if y[0] in HumanGenes:
        F_Dmelanogaster_homolog_genes.append(y)
    else:
        continue

for y in Agambiae_homolog_genes:
    if y[0] in HumanGenes:
        F_Agambiae_homolog_genes.append(y)
    else:
        continue

for y in Celegans_homolog_genes:
    if y[0] in HumanGenes:
        F_Celegans_homolog_genes.append(y)
    else:
        continue

for y in Scerevisiae_homolog_genes:
    if y[0] in HumanGenes:
        F_Scerevisiae_homolog_genes.append(y)
    else:
        continue

for y in Klactis_homolog_genes:
    if y[0] in HumanGenes:
        F_Klactis_homolog_genes.append(y)
    else:
        continue

for y in Egossypii_homolog_genes:
    if y[0] in HumanGenes:
        F_Egossypii_homolog_genes.append(y)
    else:
        continue

for y in Spombe_homolog_genes:
    if y[0] in HumanGenes:
        F_Spombe_homolog_genes.append(y)
    else:
        continue

for y in Moryzae_homolog_genes:
    if y[0] in HumanGenes:
        F_Moryzae_homolog_genes.append(y)
    else:
        continue

for y in Ncrassa_homolog_genes:
    if y[0] in HumanGenes:
        F_Ncrassa_homolog_genes.append(y)
    else:
        continue

for y in Athaliana_homolog_genes:
    if y[0] in HumanGenes:
        F_Athaliana_homolog_genes.append(y)
    else:
        continue

for y in Osativa_homolog_genes:
    if y[0] in HumanGenes:
        F_Osativa_homolog_genes.append(y)
    else:
        continue

print('Homolog genes filtered for human genes with AUG dORFs, updated lists generated')

#Convert homolog gene lists to df then save as .csv file
Hsapiens_homolog_genes_DF = pd.DataFrame({'Hsapiens_homolog_genes':F_Hsapiens_homolog_genes})
Hsapiens_homolog_genes_DF.to_csv('Hsapiens_homolog_genes.csv')
Ptroglodytes_homolog_genes_DF = pd.DataFrame({'Ptroglodytes_homolog_genes':F_Ptroglodytes_homolog_genes})
Ptroglodytes_homolog_genes_DF.to_csv('Ptroglodytes_homolog_genes.csv')
Mmulatta_homolog_genes_DF = pd.DataFrame({'Mmulatta_homolog_genes':F_Mmulatta_homolog_genes})
Mmulatta_homolog_genes_DF.to_csv('Mmulatta_homolog_genes.csv')
Clupus_homolog_genes_DF = pd.DataFrame({'Clupus_homolog_genes':F_Clupus_homolog_genes})
Clupus_homolog_genes_DF.to_csv('Clupus_homolog_genes.csv')
Btaurus_homolog_genes_DF = pd.DataFrame({'Btaurus_homolog_genes':F_Btaurus_homolog_genes})
Btaurus_homolog_genes_DF.to_csv('Btaurus_homolog_genes.csv')
Mmusculus_homolog_genes_DF = pd.DataFrame({'Mmusculus_homolog_genes':F_Mmusculus_homolog_genes})
Mmusculus_homolog_genes_DF.to_csv('Mmusculus_homolog_genes.csv')
Rnorvegicus_homolog_genes_DF = pd.DataFrame({'Rnorvegicus_homolog_genes':F_Rnorvegicus_homolog_genes})
Rnorvegicus_homolog_genes_DF.to_csv('Rnorvegicus_homolog_genes.csv')
Ggallus_homolog_genes_DF = pd.DataFrame({'Ggallus_homolog_genes':F_Ggallus_homolog_genes})
Ggallus_homolog_genes_DF.to_csv('Ggallus_homolog_genes.csv')
Xtropicalis_homolog_genes_DF = pd.DataFrame({'Xtropicalis_homolog_genes':F_Xtropicalis_homolog_genes})
Xtropicalis_homolog_genes_DF.to_csv('Xtropicalis_homolog_genes.csv')
Drerio_homolog_genes_DF = pd.DataFrame({'Drerio_homolog_genes':F_Drerio_homolog_genes})
Drerio_homolog_genes_DF.to_csv('Drerio_homolog_genes.csv')
Dmelanogaster_homolog_genes_DF = pd.DataFrame({'Dmelanogaster_homolog_genes':F_Dmelanogaster_homolog_genes})
Dmelanogaster_homolog_genes_DF.to_csv('Dmelanogaster_homolog_genes.csv')
Agambiae_homolog_genes_DF = pd.DataFrame({'Agambiae_homolog_genes':F_Agambiae_homolog_genes})
Agambiae_homolog_genes_DF.to_csv('Agambiae_homolog_genes.csv')
Celegans_homolog_genes_DF = pd.DataFrame({'Celegans_homolog_genes':F_Celegans_homolog_genes})
Celegans_homolog_genes_DF.to_csv('Celegans_homolog_genes.csv')
Scerevisiae_homolog_genes_DF = pd.DataFrame({'Scerevisiae_homolog_genes':F_Scerevisiae_homolog_genes})
Scerevisiae_homolog_genes_DF.to_csv('Scerevisiae_homolog_genes.csv')
Klactis_homolog_genes_DF = pd.DataFrame({'Klactis_homolog_genes':F_Klactis_homolog_genes})
Klactis_homolog_genes_DF.to_csv('Klactis_homolog_genes.csv')
Egossypii_homolog_genes_DF = pd.DataFrame({'Egossypii_homolog_genes':F_Egossypii_homolog_genes})
Egossypii_homolog_genes_DF.to_csv('Egossypii_homolog_genes.csv')
Spombe_homolog_genes_DF = pd.DataFrame({'Spombe_homolog_genes':F_Spombe_homolog_genes})
Spombe_homolog_genes_DF.to_csv('Spombe_homolog_genes.csv')
Moryzae_homolog_genes_DF = pd.DataFrame({'Moryzae_homolog_genes':F_Moryzae_homolog_genes})
Moryzae_homolog_genes_DF.to_csv('Moryzae_homolog_genes.csv')
Ncrassa_homolog_genes_DF = pd.DataFrame({'Ncrassa_homolog_genes':F_Ncrassa_homolog_genes})
Ncrassa_homolog_genes_DF.to_csv('Ncrassa_homolog_genes.csv')
Athaliana_homolog_genes_DF = pd.DataFrame({'Athaliana_homolog_genes':F_Athaliana_homolog_genes})
Athaliana_homolog_genes_DF.to_csv('Athaliana_homolog_genes.csv')
Osativa_homolog_genes_DF = pd.DataFrame({'Osativa_homolog_genes':F_Osativa_homolog_genes})
Osativa_homolog_genes_DF.to_csv('Osativa_homolog_genes.csv')

print('Homolog gene lists for species saved as .csv files')
