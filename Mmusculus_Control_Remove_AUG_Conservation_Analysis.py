#Python Script to updated AUG dORF conservation results by running the same analysis of AUG dORF similarity with the AUG start codons removed from the alignment
#imported required modules
import pandas as pd
import os
import re
import warnings
from pandas.errors import SettingWithCopyWarning

#prevent warning showing up every time for the loops
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

#load emboss module
os.system("module load emboss")
print('emboss module loaded successfully')

#load dataframes from .csv files and reset indexes - use the results dataframe for dORF comparisons
ResultsDF = pd.read_csv("Human_and_Mmusculus_ResultsDF_water_length_restriction_and_human_dORF_ID_duplicates_removed.csv", index_col=0)
ResultsDF = ResultsDF.reset_index(drop=True)

print('Results dataframe imported')
print('Number of lines to run analysis on:', len(ResultsDF))

#Create results dataframe with list of columns as column headings
columns = ['Human Gene', 'Human Version', 'Human 3UTR Length', 'Human 3UTRseq', 'Homolog Species', 'Homolog Gene', 'Homolog Version', 'Homolog 3UTR Length', 'Homolog 3UTRseq', 'Homolog 3UTR Length Difference', '3UTR Water Length', '3UTR Water Score', '3UTR Water Identity', '3UTR Water Similarity', '3UTR Water Gaps', 'Human dORF ID', 'Human dORFstart', 'Human dORFend', 'Human dORF Length', 'Human dORFseq', 'Homolog dORFstart', 'Homolog dORFend', 'Homolog dORF Length', 'Homolog dORFseq', 'Homolog dORF Length Difference', 'Homolog dORF Start Difference', 'Homolog dORF End Difference', 'dORF Water Length', 'dORF Water Score', 'dORF Water Identity', 'dORF Water Similarity', 'dORF Water Gaps', 'dORF vs 3UTR Water Identity Difference', 'dORF vs 3UTR Water Similarity Difference', 'dORF vs 3UTR Water Gaps Difference', 'No AUG dORF Water Length', 'No AUG dORF Water Score', 'No AUG dORF Water Identity', 'No AUG dORF Water Similarity', 'No AUG dORF Water Gaps', 'No AUG dORF vs 3UTR Water Identity Difference', 'No AUG dORF vs 3UTR Water Similarity Difference', 'No AUG dORF vs 3UTR Water Gaps Difference','No AUG dORF vs dORF Water Identity Difference', 'No AUG dORF vs dORF Water Similarity Difference','No AUG dORF vs dORF Water Gaps Difference']
ResultsDF_Ctrl = pd.DataFrame(columns=columns)

#Save Results DF as empty .csv with headings for columns
ResultsDF_Ctrl.to_csv('Human_and_Mmusculus_ResultsDF_No_AUG_dORF_water_included.csv')
print('Empty Results Dataframe .csv created')

#Create temporary subsets of results dataframes for each line in a loop
for i in range(len(ResultsDF)):
    #Create empty temp results dataframe
    tempResultsDF = pd.DataFrame(columns=columns)
    #Add row of results dataframe to temporary dataframe results
    tempResultsDF = pd.concat([tempResultsDF, ResultsDF.loc[[i]]])
    tempResultsDF = tempResultsDF.reset_index(drop=True)
    #Take human and homolog dORF sequences and take from position 4 so do not consider the AUG at the start
    HumanNoAUGdORF = tempResultsDF['Human dORFseq'][0][3:]
    HomologNoAUGdORF = tempResultsDF['Homolog dORFseq'][0][3:]
    #Write out No AUG dORF human and homolog sequences into text files to be used with water tool in command line
    Human_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'MmusHuman_temp_seq.txt'), 'w')
    Human_temp_seq.write(HumanNoAUGdORF)
    Human_temp_seq.close()
    Homolog_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'MmusHomolog_temp_seq.txt'), 'w')
    Homolog_temp_seq.write(HomologNoAUGdORF)
    Homolog_temp_seq.close()

    #run water alignment with 2 temp 3UTR sequence files
    os.system("water MmusHuman_temp_seq.txt MmusHomolog_temp_seq.txt -gapopen 10.0 -gapextend 0.5 -outfile MmusTemp_Align.water >/dev/null 2>&1")

    #remove temporary sequence files - no longer needed - save space
    os.system("rm MmusHuman_temp_seq.txt")
    os.system("rm MmusHomolog_temp_seq.txt")

    #open water alignment file from path as string in python script
    with open('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis/MmusTemp_Align.water') as water_file:
        Temp_Align = water_file.read()
    
    #Take water alignment file string line by line and pull out lines matching the results wanted    
    for lines in Temp_Align.split("\n"):
        if "Length" in lines:
            Temp_Length = lines.strip()
    for lines in Temp_Align.split("\n"):
        if "Identity" in lines:
            Temp_Identity = lines.strip()
    for lines in Temp_Align.split("\n"):
        if "Similarity" in lines:
            Temp_Similarity = lines.strip()
    for lines in Temp_Align.split("\n"):
        if "Gaps" in lines:
            Temp_Gaps = lines.strip()
    for lines in Temp_Align.split("\n"):
        if "Score" in lines:
            Temp_Score = lines.strip()

    #use regular expressions to gather water alignment values from selected lines of alignment results file
    Temp_Length = re.search("[0-9]{1,}", Temp_Length)
    Temp_Length = Temp_Length.group(0)
    Temp_Length = int(Temp_Length)
    Temp_Score = re.search("[0-9]{1,}\.[0-9]", Temp_Score)
    Temp_Score = Temp_Score.group(0)
    Temp_Score = float(Temp_Score)
    Temp_Identity = re.search("[0-9]{1,}\.[0-9]", Temp_Identity)
    Temp_Identity = Temp_Identity.group(0)
    Temp_Identity = float(Temp_Identity)
    Temp_Similarity = re.search("[0-9]{1,}\.[0-9]", Temp_Similarity)
    Temp_Similarity = Temp_Similarity.group(0)
    Temp_Similarity = float(Temp_Similarity)
    Temp_Gaps = re.search("[0-9]{1,}\.[0-9]", Temp_Gaps)
    Temp_Gaps = Temp_Gaps.group(0)
    Temp_Gaps = float(Temp_Gaps)
    
    #remove water alignment file from current directory using command line
    os.system("rm MmusTemp_Align.water")
    
    #if water alignment length is more or equal to the human or homolog temp sequence length then add results to dataframe if not not all of one of the sequences is aligned
    if Temp_Length >= len(HumanNoAUGdORF) or Temp_Length >= len(HomologNoAUGdORF):
    #add No AUG dORF water alignment results to temp df
        tempResultsDF['No AUG dORF Water Length'][0] = Temp_Length
        tempResultsDF['No AUG dORF Water Score'][0] = Temp_Score
        tempResultsDF['No AUG dORF Water Identity'][0] = Temp_Identity
        tempResultsDF['No AUG dORF Water Similarity'][0] = Temp_Similarity
        tempResultsDF['No AUG dORF Water Gaps'][0] = Temp_Gaps
        #Compare water alignments results for No AUG dORF compared to the 3' UTR between human and homolog
        threeUTRidentity = tempResultsDF['3UTR Water Identity'][0]
        tempIdentityDiff = Temp_Identity - threeUTRidentity
        tempIdentityDiff = format(tempIdentityDiff, '.2f')
        tempResultsDF['No AUG dORF vs 3UTR Water Identity Difference'] = [tempIdentityDiff]
        threeUTRsimilarity = tempResultsDF['3UTR Water Similarity'][0]
        tempSimilarityDiff = Temp_Similarity - threeUTRsimilarity
        tempSimilarityDiff = format(tempSimilarityDiff, '.2f')
        tempResultsDF['No AUG dORF vs 3UTR Water Similarity Difference'] = [tempSimilarityDiff]
        threeUTRgaps = tempResultsDF['3UTR Water Gaps'][0]
        tempGapsDiff = Temp_Gaps - threeUTRgaps
        tempGapsDiff = format(tempGapsDiff, '.2f')
        tempResultsDF['No AUG dORF vs 3UTR Water Gaps Difference'] = [tempGapsDiff]
        #compare water alignments with and without AUG included in dORFs
        dORFidentity = tempResultsDF['dORF Water Identity'][0]
        dORFtempIdentityDiff = Temp_Identity - dORFidentity
        dORFtempIdentityDiff = format(dORFtempIdentityDiff, '.2f')
        tempResultsDF['No AUG dORF vs dORF Water Identity Difference'] = [dORFtempIdentityDiff]
        dORFsimilarity = tempResultsDF['dORF Water Similarity'][0]
        dORFtempSimilarityDiff = Temp_Similarity - dORFsimilarity
        dORFtempSimilarityDiff = format(dORFtempSimilarityDiff, '.2f')
        tempResultsDF['No AUG dORF vs dORF Water Similarity Difference'] = [dORFtempSimilarityDiff]
        dORFgaps = tempResultsDF['dORF Water Gaps'][0]
        dORFtempGapsDiff = Temp_Gaps - dORFgaps
        dORFtempGapsDiff = format(dORFtempGapsDiff, '.2f')
        tempResultsDF['No AUG dORF vs dORF Water Gaps Difference'] = [dORFtempGapsDiff]
        
        #write results of temporary results dataframe into created.csv results file
        tempResultsDF.to_csv('Human_and_Mmusculus_ResultsDF_No_AUG_dORF_water_included.csv', header=None, mode='a')
        
    else:
        continue
        
    #print rolling count of number of comparisons completed
    print('\r\033[K', end='')
    print("Number of results dataframe lines analysed", i, end='\r')
    
#print line to terminal to confirm analysis completed and dataframe saved
print('\n','Analysis Complete - Dataframe saved as csv')

