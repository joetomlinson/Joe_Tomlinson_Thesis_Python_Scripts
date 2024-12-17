#Script to add dORF ID and remove duplicate 3UTRs from human dORF containing 3UTR dataframe
#imported required modules
import pandas as pd
import warnings
from pandas.errors import SettingWithCopyWarning

#prevent warining showing up every time for the loops
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

#load dataframes from .csv files and reset indexes
Human = pd.read_csv("Human_RefSeq_dORF_3UTRs.csv", index_col=0)
Human = Human.reset_index(drop=True)

print('Human dORF contatining 3UTRs imported')

#remove duplicate rows from dataframe - keep first occurence of duplicate and reset index
Human_rd = Human.drop_duplicates(subset=['Gene Name', 'Three Prime Length', 'Three Prime Sequence', 'Number of dORFs', 'dORF Locations', 'dORF Lengths', 'dORF Sequences'], keep='first', inplace=False, ignore_index=False)
Human_rd = Human_rd.reset_index(drop=True)

#remove unnamed column with previous index from.csv to df in previous script
Human_rd = Human_rd.drop('Unnamed: 0', axis=1)

print('Duplicate 3UTR sequences removed from human dataframe')

#reformat uploaded .csv dataframe - return original lists to lists from strings.
#loop through whole human dataframe
for i in range(len(Human_rd)):
    #if only one dORF treated slightly differently but convert dORF location string to list of integers
    if Human_rd['Number of dORFs'][i] == 1:
        Human_rd['dORF Locations'][i] = Human_rd['dORF Locations'][i].strip('[')
        Human_rd['dORF Locations'][i] = Human_rd['dORF Locations'][i].strip(']')
        Human_rd['dORF Locations'][i] = Human_rd['dORF Locations'][i].split(', ')
        for n in range(len(Human_rd['dORF Locations'][i])):
            Human_rd['dORF Locations'][i][n] = int(Human_rd['dORF Locations'][i][n])
    #if more than one dORF convert dORF locations string to list of list of integers
    else:
       Human_rd['dORF Locations'][i] = Human_rd['dORF Locations'][i].split('], ')
       for x in range(len(Human_rd['dORF Locations'][i])):
           Human_rd['dORF Locations'][i][x] = Human_rd['dORF Locations'][i][x].strip('[')
           Human_rd['dORF Locations'][i][x] = Human_rd['dORF Locations'][i][x].strip(']')
           Human_rd['dORF Locations'][i][x] = Human_rd['dORF Locations'][i][x].split(', ')
           for n in range(len(Human_rd['dORF Locations'][i][x])):
               Human_rd['dORF Locations'][i][x][n] = int(Human_rd['dORF Locations'][i][x][n])
    #Convert dORF lengths string into list of int lengths of dORFs
    Human_rd['dORF Lengths'][i] = Human_rd['dORF Lengths'][i].strip('[')
    Human_rd['dORF Lengths'][i] = Human_rd['dORF Lengths'][i].strip(']')
    Human_rd['dORF Lengths'][i] = Human_rd['dORF Lengths'][i].split(", ")
    for y in range(len(Human_rd['dORF Lengths'][i])):
        Human_rd['dORF Lengths'][i][y] = int(Human_rd['dORF Lengths'][i][y])
    #Convert dORF sequences str into list of dORF sequence strings
    Human_rd['dORF Sequences'][i] = Human_rd['dORF Sequences'][i].strip('[')
    Human_rd['dORF Sequences'][i] = Human_rd['dORF Sequences'][i].strip(']')
    Human_rd['dORF Sequences'][i] = Human_rd['dORF Sequences'][i].strip("'")
    Human_rd['dORF Sequences'][i] = Human_rd['dORF Sequences'][i].split("', '")
    
print('.csv upload file has been reformatted to return str variables to their original lists')

#generate dORFids first add empty column ready to hold list of dORF IDs
Human_rd.insert(12, 'Human dORF ID', '')

#Go through each line and dORF in dataset and assign dORF ID for each human dORF
for n in range(len(Human_rd)):
    Human_dORF_ID = list()
    for i in range(0, Human_rd['Number of dORFs'][n]):
        tempHuman_dORF_ID = 'dORF_' + str(i+1) + '_' + Human_rd['Accession Version'][n]
        Human_dORF_ID.append(tempHuman_dORF_ID)
    Human_rd['Human dORF ID'][n] = Human_dORF_ID
    
print('Human dORF IDs added for each dORF in dataset')

#save human dataset with ducplicate 3UTRs removed and dORF IDs added    
Human_rd.to_csv('Human_RefSeq_dORF_3UTRs_remove_duplicate_3UTR.csv')
