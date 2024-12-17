# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:30:10 2024

@author: jtoml
"""

#imported required modules
import pandas as pd
import warnings
from pandas.errors import SettingWithCopyWarning

#prevent warning showing up every time for the loops
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

#load dataframes from .csv files and reset indexes
ResultsDF = pd.read_csv("Human_and_Ggallus_ResultsDF.csv", index_col=0)
ResultsDF = ResultsDF.reset_index(drop=True)

print('Results dataframe imported')

#Create results dataframe with list of columns as column headings
columns = ['Human Gene', 'Human Version', 'Human 3UTR Length', 'Human 3UTRseq', 'Homolog Species', 'Homolog Gene', 'Homolog Version', 'Homolog 3UTR Length', 'Homolog 3UTRseq', 'Homolog 3UTR Length Difference', '3UTR Water Length', '3UTR Water Score', '3UTR Water Identity', '3UTR Water Similarity', '3UTR Water Gaps', 'Human dORF ID', 'Human dORFstart', 'Human dORFend', 'Human dORF Length', 'Human dORFseq', 'Homolog dORFstart', 'Homolog dORFend', 'Homolog dORF Length', 'Homolog dORFseq', 'Homolog dORF Length Difference', 'Homolog dORF Start Difference', 'Homolog dORF End Difference', 'dORF Water Length', 'dORF Water Score', 'dORF Water Identity', 'dORF Water Similarity', 'dORF Water Gaps', 'dORF vs 3UTR Water Identity Difference', 'dORF vs 3UTR Water Similarity Difference', 'dORF vs 3UTR Water Gaps Difference']
ResultsDF_water_restriction = pd.DataFrame(columns=columns)

#Save Results DF as empty .csv with headings for columns
ResultsDF_water_restriction.to_csv('Human_and_Ggallus_ResultsDF_water_length_restriction.csv')
print('Empty Results Dataframe .csv created')

#Create temporary subsets of results dataframes for each line in a loop
for i in range(len(ResultsDF)):
    #Create empty temp results dataframe
    tempResultsDF = pd.DataFrame(columns=columns)
    #Add row of results dataframe to temporary dataframe results
    tempResultsDF = pd.concat([tempResultsDF, ResultsDF.loc[[i]]])
    tempResultsDF = tempResultsDF.reset_index(drop=True)
    if tempResultsDF['Human 3UTR Length'][0] >= tempResultsDF['Homolog 3UTR Length'][0]:
        MIN3UTR = tempResultsDF['Homolog 3UTR Length'][0]
    else:
        MIN3UTR = tempResultsDF['Human 3UTR Length'][0]
    if tempResultsDF['Human dORF Length'][0] >= tempResultsDF['Homolog dORF Length'][0]:
        MINdORF = tempResultsDF['Homolog dORF Length'][0]
    else:
        MINdORF = tempResultsDF['Human dORF Length'][0]
    if tempResultsDF['3UTR Water Length'][0] >= MIN3UTR and tempResultsDF['dORF Water Length'][0] >= MINdORF:
        #write results of temporary results dataframe into created.csv results file
        tempResultsDF.to_csv('Human_and_Ggallus_ResultsDF_water_length_restriction.csv', header=None, mode='a')
    else:
        continue
    
    #print rolling count of number of comparisons completed
    print('\r\033[K', end='')
    print("Number of results dataframe lines analysed", i, end='\r')
    
#print line to terminal to confirm analysis completed and dataframe saved
print('\n','Analysis Complete - Dataframe saved as csv')

#load dataframes from .csv files and reset indexes
ResultsDF_no_dup = pd.read_csv("Human_and_Ggallus_ResultsDF_water_length_restriction.csv", index_col=0)
ResultsDF_no_dup = ResultsDF_no_dup.reset_index(drop=True)

print('Results dataframe imported')

#Sort dataframe in dec order in terms of dORF similarity %
ResultsDF_no_dup.sort_values(by=['dORF Water Similarity'], ascending=False)

print('Results dataframe sorted descending for dORF similarity')

#remove duplicate rows from dataframe - keep first occurence of duplicate and reset index
ResultsDF_no_dup = ResultsDF_no_dup.drop_duplicates(subset=['Human dORF ID'], keep='first', inplace=False, ignore_index=False)
ResultsDF_no_dup = ResultsDF_no_dup.reset_index(drop=True)

print('Duplicate Human dORF IDs removed from results')

#save results dataframe with duplicate human dORF IDs removed    
ResultsDF_no_dup.to_csv('Human_and_Ggallus_ResultsDF_water_length_restriction_and_human_dORF_ID_duplicates_removed.csv')

print('Updated Results dataframe saved as .csv')

