#Python Script to take Human RefSeq Accessions and gather all 3' UTR sequences from NCBI nucleotide database then split into dataframes with and without AUG dORFs in these sequences, also include details about dORFs in dataframe

import pandas as pd #Import pandas
import re #import regular expressions package
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #prevent FutureWarning note printed about df.append for each query

df = pd.read_csv("Human_Accession_Versions.csv") #convert human accession csv file in current directory to dataframe

Human_ID_List = df['Accession Version'].to_list() #Convert all column 1 values into list of human accession versions with refseq
Human_ID_List_1 = Human_ID_List[0:50000]
Human_ID_List_2 = Human_ID_List[50000:100000]
Human_ID_List_3 = Human_ID_List[100000:150000]
Human_ID_List_4 = Human_ID_List[150000:200000]

print("Human_ID_List Generated") #Show this at command line

columns = ['Accession Version', 'Accession', 'Gene Name', 'Gene Description', 'Transcript Annotation', 'Species', 'Transcript Sequence', 'Three Prime Start', 'Three Prime End', 'Three Prime Length', 'Three Prime Sequence'] #Create list of column headings

Human_CDS_end_not_integer = list() #Create empty list to store UIDs of transcripts which have no CDS end integer
Human_no_3primeUTR_collected = list() #Create empty list to store UIDs of transcripts which have no 3' UTR

#Carry out script for Part 1 of human IDs
df_RefSeq_Three_prime = pd.DataFrame(columns=columns) #Create empty Table with headings from list ready for data for each accession searched

from Bio import Entrez #Use Entrez toolpack from biopython

Entrez.email = "j.tomlinson@mail.com" #Provide NCBI Entrez with email
Entrez.api_key = "a3427acfca962bf95e67b352748ce9154508" #Provide NCBI with API increase requests per second 3-10

#Loop to go through all 199205 human IDs but splitting them into 200 at a time to feed to Entrez.efetch to avoid overloading it and speed things up
for n in range(0, len(Human_ID_List_4), 200):
    while True: #while true, try, except and break brings in a forcing feature where the loop is forces to keep rerunning when it fails until it runns successfully
        try:
            TEMP_Human_IDs = Human_ID_List_4[n:n+200] #Create temporary list to store the 200 at a time human IDs to use in next line
            Entrez_handle = Entrez.efetch(db="nucleotide", id=TEMP_Human_IDs, rettype="gb", retmode="xml") #Search Entrez NCBI nucleotide database for selected 200 accessions and generate handle to access genbank data for this accession
            Entrez_headings = Entrez.read(Entrez_handle) #Access full genebank record and read - produces list of 200 records - for each accession there is the genbank record
        except:
            continue #part of loop to rerun the script if the server has failed 
        break

    #Create for loop go through each of the 200 genbank records in list and gather data for each in the overall table
    for i in range(len(Entrez_headings)):
        TEMP_df_RefSeq_Three_prime = pd.DataFrame(columns=columns) #Create empty temporary table same as previous to be used in for loop to be refilled with each accession as it goes
                
        TEMP_Entrez_headings = Entrez_headings[i] #Access full genebank record and read - one record (ID) at a time following through loop
                
        Gene_name = str(TEMP_Entrez_headings['GBSeq_definition']) #Create string variable with all data from gene definition
        Gene_name = Gene_name.split("(") #split spring with ( and replace gene name variable
        Gene_name = str(Gene_name[1]) #replace gene name variable with second part of split after (
        Gene_name = Gene_name.split(")") #split string and rename variable with all before )
        Gene_name = Gene_name[0] #gene name variable with gene name taken from () in original data
        TEMP_df_RefSeq_Three_prime['Gene Name'] = [Gene_name] #Add Gene_Name variable to temp table
                
        Gene = str(TEMP_Entrez_headings['GBSeq_definition']) #Create variable with description of gene
        TEMP_df_RefSeq_Three_prime['Gene Description'] = [Gene] #Add Gene variable to temp table
        Accession = str(TEMP_Entrez_headings['GBSeq_locus']) #Create variable with accession code
        TEMP_df_RefSeq_Three_prime['Accession'] = [Accession] #Add Accession variable to temp table
        Accession_version = str(TEMP_Entrez_headings['GBSeq_accession-version']) #Create variable with accession version ie specific to transcript version
        TEMP_df_RefSeq_Three_prime['Accession Version'] = [Accession_version] #Add Accession_version variable to temp table
        Species = str(TEMP_Entrez_headings['GBSeq_organism']) #Create variable with species name from genbank entry
        TEMP_df_RefSeq_Three_prime['Species'] = [Species] #Add Species variable to temp table
        Transcript_Annotation = str(TEMP_Entrez_headings['GBSeq_keywords']) #Create variable with keyword annotation from genbank entry
        TEMP_df_RefSeq_Three_prime['Transcript Annotation'] = [Transcript_Annotation] #Add Transcript_Annotation variable to temp table
        Transcript_sequence = str(TEMP_Entrez_headings['GBSeq_sequence']) #Create variable with transcript sequence from genbank entry
        TEMP_df_RefSeq_Three_prime['Transcript Sequence'] = [Transcript_sequence] #Add transcript sequence variable to temp table
                
        Entrez_features = str(TEMP_Entrez_headings['GBSeq_feature-table']) #Create string variable with all genbank data related to features table
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
            Human_CDS_end_not_integer.append([[Gene_name],[Accession_version]]) #store the errors - keep list in list of human gene, UID with no CDS integer
            continue
            
        Three_prime_start = CDS_end + 1 #create variable with 3UTR start which is 1 base further than CDS end
        TEMP_df_RefSeq_Three_prime['Three Prime Start'] = [Three_prime_start] #Add 3UTR start variable to temp table
        
        Three_prime_end = int(TEMP_Entrez_headings['GBSeq_length']) #create integer variable as end of 3UTR same as total length
        TEMP_df_RefSeq_Three_prime['Three Prime End'] = [Three_prime_end] #Add 3UTR end variable to temp table
        
        Three_prime_sequence = str(Transcript_sequence[(Three_prime_start-1):]) #Create variable with 3UTR sequence by taking 3UTR start site adjusting for locations starting at 0 by doing -1 to the end of the transcript
        pattern = r'[^\.AaCcGgTtUuXxNn]' #pattern to search for any character other then ACGTUXN regardless of case
        if re.search(pattern, Three_prime_sequence):
            Human_no_3primeUTR_collected.append([[Gene_name],[Accession_version]]) #store the errors - keep list in list of human gene, UID with no 3' UTR
            continue #If there is no 3'UTR sequence gathered move onto next i value dont add to larger 3' UTR table
        else:
            TEMP_df_RefSeq_Three_prime['Three Prime Sequence'] = [Three_prime_sequence] #Add 3UTR sequence variable to temp table
        
        Three_prime_length = len(Three_prime_sequence) #store length of 3UTR sequence as variable
        if Three_prime_length == 0:
            Human_no_3primeUTR_collected.append([[Gene_name],[Accession_version]]) #store the errors - keep list in list of human gene, UID with no 3' UTR
            continue #If there is no 3'UTR sequence gathered move onto next i value dont add to larger 3' UTR table
        TEMP_df_RefSeq_Three_prime['Three Prime Length'] = [Three_prime_length] #Add 3UTR length variable to temp table
        
        #Add data from temp table to larger overall table each for loop and accession in new rows for each accession
        df_RefSeq_Three_prime = pd.concat([df_RefSeq_Three_prime, TEMP_df_RefSeq_Three_prime])
        df_RefSeq_Three_prime = df_RefSeq_Three_prime.reset_index(drop=True)
        
        print('\r\033[K', end='')
        print("Human Accessions Completed -", n, end='\r') #Print number of accessions run on same line as rolling count - clear then print 
print("Human Accessions Analysis Complete")

print("Human transcripts with CDS end not integer -", len(Human_CDS_end_not_integer))
print("Human transcripts with no 3' UTR sequence collected -", len(Human_no_3primeUTR_collected))

del (Accession, Accession_version, CDS_end, CDS_feature, CDS_location, columns, df, Entrez_features, Entrez_features_list, Entrez_handle, Entrez_headings, Gene, Gene_name, Human_ID_List, i, n, Species, TEMP_df_RefSeq_Three_prime, TEMP_Entrez_headings, TEMP_Human_IDs, Three_prime_end, Three_prime_length, Three_prime_sequence, Three_prime_start, Transcript_Annotation, Transcript_sequence) #delete no longer needed variables

###Find dORF sequences and add to dataframe###
df_RefSeq_Three_prime["Number of dORFs"] = "" #Add new empty column to dataframe ready for number of dORFs
df_RefSeq_Three_prime["dORF Locations"] = "" #Add new empty column to dataframe ready for dORF locations
df_RefSeq_Three_prime["dORF Lengths"] = "" #Add new empty column to dataframe ready for dORF sequence lengths
df_RefSeq_Three_prime["dORF Sequences"] = "" #Add new empty column to dataframe ready for dORF sequences

import re #import regular expressions package

regex = re.compile(r'(?=(atg(?:...)*?(tag|tga|taa)))', re.IGNORECASE) #set regular expression to find all ATG ORFs - changed to lowercase - added code so that the regular expression is case insensitive

for i in range(len(df_RefSeq_Three_prime)): #loop to go through each 3' UTR sequence
     Temp_3UTR_seq = df_RefSeq_Three_prime['Three Prime Sequence'][i] #Create temporary variable with each 3' UTR sequence 
     Temp_rawAllORFs = regex.findall(Temp_3UTR_seq) #find all str matching to regular expression in Seq variable
     Temp_dORFs = list() #creat empty list to add dORF sequences to later
     
     for n in range(len(Temp_rawAllORFs)): #create loop to go through each ORF found
         tempORF = Temp_rawAllORFs[n][0] #only keep ORF sequence not stop codon on its own
         if len(tempORF)>29 and len(tempORF)<301: #only add dORF sequence to list if its length is between 30-300 bases
             Temp_dORFs.append(tempORF)
         else:
             continue #move onto next ORF entry in list if above if statement not met
             
     df_RefSeq_Three_prime['dORF Sequences'][i] = Temp_dORFs #add dORF sequences to dataframe
     df_RefSeq_Three_prime['Number of dORFs'][i] = len(Temp_dORFs) #add number of dORF sequences to dataframe
     
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

     df_RefSeq_Three_prime['dORF Locations'][i] = Temp_dORFlocations #add dORF locations to dataframe
     df_RefSeq_Three_prime['dORF Lengths'][i] = Temp_dORFlengths #add dORF lengths to dataframe
     
     print('\r\033[K', end='')
     print("3' UTR sequences searched for dORF -", i, end='\r') #Print number of accessions run on same line as rolling count - clear then print 
del(i, n, x, regex, Temp_3UTR_seq, Temp_rawAllORFs, Temp_dORFs, tempORF, Temp_dORFlocations, Temp_dORFlengths,TempdORFstart, TempdORFend, TempdORFlength) #remove any no longer wanted variables
print("Human 3' UTR sequences search for dORFs - Analysis Complete")
#save dataframe as csv in current directory
df_RefSeq_Three_prime.to_csv('df_RefSeq_Three_prime_4.csv')

#Split above dataframe into human 3' UTRs with and without dORFs
#Human_RefSeq_dORF_3UTRs = df_RefSeq_Three_prime[df_RefSeq_Three_prime['Number of dORFs'] > 0]
#Human_RefSeq_NO_dORF_3UTRs = df_RefSeq_Three_prime[df_RefSeq_Three_prime['Number of dORFs'] == 0]

#save dataframe as csv in current directory
#Human_RefSeq_dORF_3UTRs.to_csv('Human_RefSeq_dORF_3UTRs.csv')
#Human_RefSeq_NO_dORF_3UTRs.to_csv('Human_RefSeq_NO_dORF_3UTRs.csv')

#Convert lists of missing values to dataframes and then save as csv files in current directory
Human_CDS_end_not_integer_DF = pd.DataFrame({'Human_CDS_end_not_integer':Human_CDS_end_not_integer})
Human_no_3primeUTR_collected_DF = pd.DataFrame({'Human_no_3primeUTR_collected':Human_no_3primeUTR_collected})
Human_CDS_end_not_integer_DF.to_csv('Human_CDS_end_not_integer_4.csv')
Human_no_3primeUTR_collected_DF.to_csv('Human_no_3primeUTR_collected_4.csv')

print("Human Dataframes with and without dORF 3UTRs saved as csv")
