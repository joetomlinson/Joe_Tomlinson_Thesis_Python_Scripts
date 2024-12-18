# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 10:07:49 2024

@author: jtoml
"""

#### Example script for Ncrassa to collect nucleotide UIDs and 3UTR and dORF sequences for homolog genes would need to run for each species ####

#Import pandas and prevent FutureWarning note printed about df.append for each query
import pandas as pd 
import re
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#Use Entrez toolpack from biopython and provide NCBI Entrez with email and API - increase requests per second 3-10
from Bio import Entrez 
Entrez.email = "j.tomlinson@mail.com"
Entrez.api_key = "a3427acfca962bf95e67b352748ce9154508"

#convert homolog gene csv file in current directory to dataframe and return previous list values back to list from string due to .csv file saving
df_list = pd.read_csv("Ncrassa_homolog_genes.csv", index_col=0) 
for i in range(len(df_list)):
    df_list['Ncrassa_homolog_genes'][i] = df_list['Ncrassa_homolog_genes'][i].strip("[")
    df_list['Ncrassa_homolog_genes'][i] = df_list['Ncrassa_homolog_genes'][i].strip("]")
    df_list['Ncrassa_homolog_genes'][i] = df_list['Ncrassa_homolog_genes'][i].strip("'")
    df_list['Ncrassa_homolog_genes'][i] = df_list['Ncrassa_homolog_genes'][i].split("', '")

print('Homolog genes loaded and reformated, number of homolog genes:', len(df_list))

#Create empty list to store UID and nucIDs collected by the next loop
Ncrassa_UID_or_nucID = list()

#Loop through the homolog dataframe for each homolog gene collected collect the nucID, homolog gene and human gene, use homolog gene to complete search term for NCBI nucleotide database
for i in range(len(df_list)):
    Ncrassa_nucleotide_id = df_list['Ncrassa_homolog_genes'][i][2]
    if (Ncrassa_nucleotide_id.__contains__("M")):
        Ncrassa_gene = df_list['Ncrassa_homolog_genes'][i][1]
        Human_Gene = df_list['Ncrassa_homolog_genes'][i][0]
    else:
        continue
    
    Temp_search_term = str("refseq[filter] AND biomol_mrna[PROP] AND Neurospora Crassa[organism] AND ") + Ncrassa_gene #create search term variable with all necessary info to restrict organism, refseq and mrna

    #while true, try, except and break brings in a forcing feature where the loop is forced to keep rerunning when it fails until it runs successfully - helps with NCBI server crashes
    while True: 
        try:
            #search NCBI nucleotide database with search term to get the nucleotide UID for transcripts to use with efetch
            NcrassaNucleotideUID_eSearch = Entrez.esearch(db="nucleotide", term=Temp_search_term, retmode="xml") 
            #create variable with output from esearch in digestible format - only really the ID included
            NcrassaNucleotideUID_eSearch_Results = Entrez.read(NcrassaNucleotideUID_eSearch) 
        except:
            continue
        break
    
    #Take what's included in IdList section of esearch results and store as list variable
    Ncrassa_UID_entrez = NcrassaNucleotideUID_eSearch_Results['IdList'] 
    #only add to the nucleotide UID or list of UIDs to list if there is a UID found - not empty
    if len(Ncrassa_UID_entrez) > 0: 
        #likely to be multiple IDs from one gene search for different transcripts needs loop to go through each, and Take ID part of the dictionary variable as a string - the UID from nucleotide
        for n in range(len(Ncrassa_UID_entrez)): 
            Temp_Ncrassa_UID = str(Ncrassa_UID_entrez[n])
            #if UID is not all digits then use nucleotide id instead of UID - won't run with NCBI nucleotide
            if Temp_Ncrassa_UID.isdigit():
                Ncrassa_UID_or_nucID.append([[Human_Gene],[Temp_Ncrassa_UID],[Ncrassa_gene]])
            else:
                Ncrassa_UID_or_nucID.append([[Human_Gene],[Ncrassa_nucleotide_id],[Ncrassa_gene]]) 
    #If there is no UID generated from search with gene name use the nucleotide ID from homologene database output instead add this to list
    else:
        Ncrassa_UID_or_nucID.append([[Human_Gene],[Ncrassa_nucleotide_id],[Ncrassa_gene]]) 
        continue
    
    #Print number of genes run on same line as rolling count - clear then print
    print('\r\033[K', end='')
    print("Ncrassa Homolog genes run -", i, end='\r')  

#Take list of Ncrassa nucleotide UIDs and store as .csv by converting to dataframe
Ncrassa_Homolog_UID_or_nucID_List_DF = pd.DataFrame({'Ncrassa_Homolog_UID_or_nucID_List':Ncrassa_UID_or_nucID})
Ncrassa_Homolog_UID_or_nucID_List_DF.to_csv('Ncrassa_Homolog_UID_or_nucID_List.csv')

print('Homolog UID or nucID collection complete, file saved, number of IDs:', len(Ncrassa_UID_or_nucID))

#create list of column headings for dataframe
columns = ['Human Gene Name', 'Accession Version', 'Accession', 'Gene Name', 'Gene Description', 'Transcript Annotation', 'Species', 'Transcript Sequence', 'Three Prime Start', 'Three Prime End', 'Three Prime Length', 'Three Prime Sequence'] #Create list of column headings

#create homolog dataframe to collect 3UTR info with column headings defined
Ncrassa_df_RefSeq_Three_prime = pd.DataFrame(columns=columns) #Create empty Table with headings from list ready for data for each accession searched

#Create empty list to store IDs not being entered into dataframe as gene matched is wrong, of transcripts which have no CDS end integer, and of transcripts which have no 3' UTR
Ncrassa_nucleotide_UID_not_matching_search_gene = list() 
Ncrassa_CDS_end_not_integer = list() 
Ncrassa_no_3primeUTR_collected = list()

#Loop through each ID in homolog list generated above with each entries gene, id and human gene assigned as variables below
for n in range(len(Ncrassa_UID_or_nucID)):
    Ncrassa_id = str(Ncrassa_UID_or_nucID[n][1])
    Ncrassa_id = Ncrassa_id.strip("[")
    Ncrassa_id = Ncrassa_id.strip("]")
    Ncrassa_id = Ncrassa_id.strip("'")
    Ncrassa_gene = str(Ncrassa_UID_or_nucID[n][2])
    Ncrassa_gene = Ncrassa_gene.strip("[")
    Ncrassa_gene = Ncrassa_gene.strip("]")
    Ncrassa_gene = Ncrassa_gene.strip("'")
    Human_Gene = str(Ncrassa_UID_or_nucID[n][0])
    Human_Gene = Human_Gene.strip("[")
    Human_Gene = Human_Gene.strip("]")
    Human_Gene = Human_Gene.strip("'")
    
    #while true, try, except and break brings in a forcing feature where the loop is forced to keep rerunning when it fails until it runs successfully
    while True: 
        try:
            #Search Entrez NCBI nucleotide database for id and generate handle to access genbank data for this accession and create variable with output from esearch in digestible format
            Ncrassa_Entrez_handle = Entrez.efetch(db="nucleotide", id=Ncrassa_id, rettype="gb", retmode="xml") 
            Ncrassa_Entrez_headings = Entrez.read(Ncrassa_Entrez_handle)
        except:
            continue
        break
    
    #loop through list of genbank data gathered by search
    for i in range(len(Ncrassa_Entrez_headings)):
        #Create empty temporary table same as previous to be used in for loop to be refilled with each accession as it goes
        TEMP_Ncrassa_df_RefSeq_Three_prime = pd.DataFrame(columns=columns) 
        #create variable storing each list value
        TEMP_Ncrassa_Entrez_headings = Ncrassa_Entrez_headings[i]
        
        TEMP_Ncrassa_df_RefSeq_Three_prime['Human Gene Name'] = [Human_Gene] #add human gene name to column of dataframe
        Gene_name = str(TEMP_Ncrassa_Entrez_headings['GBSeq_definition']) #Create string variable with all data from gene definition
        if 'M' in Ncrassa_id:
            TEMP_Ncrassa_df_RefSeq_Three_prime['Gene Name'] = [Ncrassa_gene]
        else:
            if (Gene_name.__contains__(Ncrassa_gene)): #really important step if the nucleotide UID does not return a result for the gene name used to search for the UID do not add to dataframe
                TEMP_Ncrassa_df_RefSeq_Three_prime['Gene Name'] = [Ncrassa_gene] #Add Gene_Name variable to temp table
            else:
                Ncrassa_nucleotide_UID_not_matching_search_gene.append([[Gene_name],[Ncrassa_gene],[Ncrassa_id]]) #store the errors - list in list store the UIDs that call the wrong gene with human gene, UID, searched homolog gene and wrongly called gene
                continue
        
        Gene = str(TEMP_Ncrassa_Entrez_headings['GBSeq_definition']) #Create variable with description of gene
        TEMP_Ncrassa_df_RefSeq_Three_prime['Gene Description'] = [Gene] #Add Gene variable to temp table
        Accession = str(TEMP_Ncrassa_Entrez_headings['GBSeq_locus']) #Create variable with accession code
        TEMP_Ncrassa_df_RefSeq_Three_prime['Accession'] = [Accession] #Add Accession variable to temp table
        Accession_version = str(TEMP_Ncrassa_Entrez_headings['GBSeq_accession-version']) #Create variable with accession version ie specific to transcript version
        TEMP_Ncrassa_df_RefSeq_Three_prime['Accession Version'] = [Accession_version] #Add Accession_version variable to temp table
        Species = str(TEMP_Ncrassa_Entrez_headings['GBSeq_organism']) #Create variable with species name from genbank entry
        TEMP_Ncrassa_df_RefSeq_Three_prime['Species'] = [Species] #Add Species variable to temp table
        Transcript_Annotation = str(TEMP_Ncrassa_Entrez_headings['GBSeq_keywords']) #Create variable with keyword annotation from genbank entry
        TEMP_Ncrassa_df_RefSeq_Three_prime['Transcript Annotation'] = [Transcript_Annotation] #Add Transcript_Annotation variable to temp table
        Transcript_sequence = str(TEMP_Ncrassa_Entrez_headings['GBSeq_sequence']) #Create variable with transcript sequence from genbank entry
        TEMP_Ncrassa_df_RefSeq_Three_prime['Transcript Sequence'] = [Transcript_sequence] #Add transcript sequence variable to temp table
                
        Entrez_features = str(TEMP_Ncrassa_Entrez_headings['GBSeq_feature-table']) #Create string variable with all genbank data related to features table
        Entrez_features_list = Entrez_features.split("'GBFeature_key': ") #Split previous string into list of features
        CDS_feature = [i for i in Entrez_features_list if i.startswith("'CDS'")] #Create variable from item in list of features that starts with CDS the CDS features section
        CDS_feature = str(CDS_feature) #convert variable into string
        CDS_feature = CDS_feature.split(" ") #split string where there are spaces into new list of same name
        CDS_location = str(CDS_feature[2]) #take 3rd entry in list which contains CDS location as new variable
        CDS_location = CDS_location.strip(",") #remove "," from string variable
        CDS_location = CDS_location.strip("'") #remove "'" from string variable allows it to be an integer later
        CDS_location = CDS_location.strip("\\") #remove "\" from string variable allows it to be an integer later
        CDS_location = CDS_location.split("..") #split string into list of CDS start and stop at .. point
        try:
            CDS_end = int(CDS_location[1]) #take second list entry as integer as variable called CDS end
        except ValueError: #If there is an error with the above line and CDS_end is not an integer then a 3'UTR sequence can't be generated so skip to next i accession and dont add this i values findings to 3'UTR table
            Ncrassa_CDS_end_not_integer.append([[Human_Gene],[Ncrassa_gene],[Ncrassa_id]]) #store the errors - keep list in list of human gene, UID and ptroglodyte gene with no CDS integer
            continue
            
        Three_prime_start = CDS_end + 1 #create variable with 3UTR start which is 1 base further than CDS end
        TEMP_Ncrassa_df_RefSeq_Three_prime['Three Prime Start'] = [Three_prime_start] #Add 3UTR start variable to temp table
        
        Three_prime_end = int(TEMP_Ncrassa_Entrez_headings['GBSeq_length']) #create integer variable as end of 3UTR same as total length
        TEMP_Ncrassa_df_RefSeq_Three_prime['Three Prime End'] = [Three_prime_end] #Add 3UTR end variable to temp table
        
        Three_prime_sequence = str(Transcript_sequence[(Three_prime_start-1):]) #Create variable with 3UTR sequence by taking 3UTR start site adjusting for locations starting at 0 by doing -1 to the end of the transcript
        pattern = r'[^\.AaCcGgTtUuXxNn]' #pattern to search for any character other then ACGTUXN regardless of case
        if re.search(pattern, Three_prime_sequence):
            Ncrassa_no_3primeUTR_collected.append([[Human_Gene],[Ncrassa_gene],[Ncrassa_nucleotide_id]]) #store the errors - keep list in list of human gene, UID and ptroglodyte gene with no 3' UTR
            continue #If there is no 3'UTR sequence gathered move onto next i value dont add to larger 3' UTR table
        else:
            TEMP_Ncrassa_df_RefSeq_Three_prime['Three Prime Sequence'] = [Three_prime_sequence] #Add 3UTR sequence variable to temp table
        
        Three_prime_length = len(Three_prime_sequence) #store length of 3UTR sequence as variable
        if Three_prime_length == 0:
            Ncrassa_no_3primeUTR_collected.append([[Human_Gene],[Ncrassa_gene],[Ncrassa_id]]) #store the errors - keep list in list of human gene, UID and ptroglodyte gene with no 3' UTR
            continue #If there is no 3'UTR sequence gathered move onto next i value dont add to larger 3' UTR table
        TEMP_Ncrassa_df_RefSeq_Three_prime['Three Prime Length'] = [Three_prime_length] #Add 3UTR length variable to temp table
                
        #Add data from temp table to larger overall table each for loop and accession in new rows for each accession
        Ncrassa_df_RefSeq_Three_prime = pd.concat([Ncrassa_df_RefSeq_Three_prime, TEMP_Ncrassa_df_RefSeq_Three_prime])
        Ncrassa_df_RefSeq_Three_prime = Ncrassa_df_RefSeq_Three_prime.reset_index(drop=True)
        
        print('\r\033[K', end='')
        print("Ncrassa UID or nucID Accessions Run -", n, end='\r') #Print number of UID accessions run per 100 on same line as rolling count - clear then print 

print("Ncrassa 3' UTR collection Analysis Complete, 3UTRs collected:", len(Ncrassa_df_RefSeq_Three_prime))
print("Ncrassa nucleotide UIDs not matching searched gene -", len(Ncrassa_nucleotide_UID_not_matching_search_gene))
print("Ncrassa nucleotide UIDs where CDS end is not integer -", len(Ncrassa_CDS_end_not_integer))
print("Ncrassa nucleotide UIDs where 3' UTR could not be collected -", len(Ncrassa_no_3primeUTR_collected))

#Convert lists of missing values to dataframes and then save as csv files in current directory
Ncrassa_nucleotide_UID_not_matching_search_gene_DF = pd.DataFrame({'Ncrassa_nucleotide_UID_not_matching_search_gene':Ncrassa_nucleotide_UID_not_matching_search_gene})
Ncrassa_CDS_end_not_integer_DF = pd.DataFrame({'Ncrassa_CDS_end_not_integer':Ncrassa_CDS_end_not_integer})
Ncrassa_no_3primeUTR_collected_DF = pd.DataFrame({'Ncrassa_no_3primeUTR_collected':Ncrassa_no_3primeUTR_collected})
Ncrassa_nucleotide_UID_not_matching_search_gene_DF.to_csv('Ncrassa_nucleotide_UID_not_matching_search_gene.csv')
Ncrassa_CDS_end_not_integer_DF.to_csv('Ncrassa_CDS_end_not_integer.csv')
Ncrassa_no_3primeUTR_collected_DF.to_csv('Ncrassa_no_3primeUTR_collected.csv')


###Find dORF sequences and add to dataframe###
Ncrassa_df_RefSeq_Three_prime["Number of dORFs"] = "" #Add new empty column to dataframe ready for number of dORFs
Ncrassa_df_RefSeq_Three_prime["dORF Locations"] = "" #Add new empty column to dataframe ready for dORF locations
Ncrassa_df_RefSeq_Three_prime["dORF Lengths"] = "" #Add new empty column to dataframe ready for dORF sequence lengths
Ncrassa_df_RefSeq_Three_prime["dORF Sequences"] = "" #Add new empty column to dataframe ready for dORF sequences

import re #import regular expressions package

regex = re.compile(r'(?=(atg(?:...)*?(tag|tga|taa)))', re.IGNORECASE) #set regular expression to find all ATG ORFs - changed to lowercase - added code so that the regular expression is case insensitive

for i in range(len(Ncrassa_df_RefSeq_Three_prime)): #loop to go through each 3' UTR sequence
     Temp_3UTR_seq = Ncrassa_df_RefSeq_Three_prime['Three Prime Sequence'][i] #Create temporary variable with each 3' UTR sequence 
     Temp_rawAllORFs = regex.findall(Temp_3UTR_seq) #find all str matching to regular expression in Seq variable
     Temp_dORFs = list() #creat empty list to add dORF sequences to later
     
     for n in range(len(Temp_rawAllORFs)): #create loop to go through each ORF found
         tempORF = Temp_rawAllORFs[n][0] #only keep ORF sequence not stop codon on its own
         if len(tempORF)>29 and len(tempORF)<301: #only add dORF sequence to list if its length is between 30-300 bases
             Temp_dORFs.append(tempORF)
         else:
             continue #move onto next ORF entry in list if above if statement not met
             
     Ncrassa_df_RefSeq_Three_prime['dORF Sequences'][i] = Temp_dORFs #add dORF sequences to dataframe
     Ncrassa_df_RefSeq_Three_prime['Number of dORFs'][i] = len(Temp_dORFs) #add number of dORF sequences to dataframe
     
     Temp_dORFlocations = list() #create empty list to store the start and stop locations of dORF sequences in Seq 
     Temp_dORFlengths = list() #create empty list to store the length of dORF sequences in Seq
     
     for x in range(len(Temp_dORFs)): #create loop to go through each dORF sequence found
         TempdORFstart = re.search(Temp_dORFs[x], Temp_3UTR_seq).start() #search for dORF sequence in Seq and report start index in string store as temporary variable
         TempdORFstart = TempdORFstart + 1 #Add 1 to dORF start site due to python index starting at 0
         TempdORFend = re.search(Temp_dORFs[x], Temp_3UTR_seq).end() #search for dORF sequence in Seq and report end index in string store as temporary variable
         TempdORFend = TempdORFend + 1 #Add 1 to dORF end site due to python index starting at 0
         Temp_dORFlocations.append([TempdORFstart, TempdORFend]) #Add dORF start and end locations as list oin list
         TempdORFlength = len(Temp_dORFs[x]) #store length of dORF sequence as temporary variable
         Temp_dORFlengths.append(TempdORFlength) #Add dORF length to list

     Ncrassa_df_RefSeq_Three_prime['dORF Locations'][i] = Temp_dORFlocations #add dORF locations to dataframe
     Ncrassa_df_RefSeq_Three_prime['dORF Lengths'][i] = Temp_dORFlengths #add dORF lengths to dataframe
     
     print('\r\033[K', end='')
     print("3' UTR sequences searched for dORF -", i, end='\r') #Print number of accessions run on same line as rolling count - clear then print 
print("Ncrassa 3' UTR sequences search for dORFs - Analysis Complete")

#save dataframe as csv in current directory
Ncrassa_df_RefSeq_Three_prime.to_csv('Ncrassa_Homologs_of_Human_RefSeq_dORF_3UTRs.csv')

print('Full Homolog dataframe saved')

#Split above dataframe into Ncrassa 3' UTRs with and without dORFs
Ncrassa_RefSeq_dORF_containing_3UTRs = Ncrassa_df_RefSeq_Three_prime[Ncrassa_df_RefSeq_Three_prime['Number of dORFs'] > 0]
Ncrassa_RefSeq_NO_dORF_3UTRs = Ncrassa_df_RefSeq_Three_prime[Ncrassa_df_RefSeq_Three_prime['Number of dORFs'] == 0]


#save dataframe as csv in current directory
Ncrassa_RefSeq_dORF_containing_3UTRs.to_csv('Ncrassa_RefSeq_dORF_containing_3UTRs.csv')
Ncrassa_RefSeq_NO_dORF_3UTRs.to_csv('Ncrassa_RefSeq_NO_dORF_3UTRs.csv')

print('Homolog dataframe split into those with and without dORFs in 3UTR and saved, entries with dORF:', len(Ncrassa_RefSeq_dORF_containing_3UTRs), '  entries without:', len(Ncrassa_RefSeq_NO_dORF_3UTRs))

#Open homolog dataframe .csv which has dORF containing 3' UTRs
Ncrassa_RefSeq_dORF_containing_3UTRs_csv = pd.read_csv("Ncrassa_RefSeq_dORF_containing_3UTRs.csv", index_col=0)
Ncrassa_RefSeq_dORF_containing_3UTRs_csv = Ncrassa_RefSeq_dORF_containing_3UTRs_csv.reset_index(drop=True)

#remove duplicate rows from dataframe - keep first occurence of duplicate and reset index and save as .csv
Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs = Ncrassa_RefSeq_dORF_containing_3UTRs_csv.drop_duplicates(subset=['Gene Name', 'Three Prime Length', 'Three Prime Sequence', 'Number of dORFs', 'dORF Locations', 'dORF Lengths', 'dORF Sequences'], keep='first', inplace=False, ignore_index=False)
Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs = Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs.reset_index(drop=True)
Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs.to_csv('Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs.csv')

print('Duplicate 3UTRs removed from homolog dataframe with dORFs, entries remaining:', len(Ncrassa_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs))
