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

#load dataframes from .csv files and reset indexes - use the dataframe with removed duplicate 3UTR sequences
ResultsDF = pd.read_csv("Human_and_Xtropicalis_ResultsDF_duplicate_3UTR_comparisons_removed.csv", index_col=0)
ResultsDF = ResultsDF.reset_index(drop=True)

print('Results dataframe imported')
print('Number of lines to run analysis on:', len(ResultsDF))

#Create results dataframe with list of columns as column headings
columns = ['Human Gene', 'Human Version', 'Human 3UTR Length', 'Human 3UTRseq', 'Homolog Species', 'Homolog Gene', 'Homolog Version', 'Homolog 3UTR Length', 'Homolog 3UTRseq', 'Homolog 3UTR Length Difference', '3UTR Water Length', '3UTR Water Score', '3UTR Water Identity', '3UTR Water Similarity', '3UTR Water Gaps', 'Human dORF ID', 'Human dORFstart', 'Human dORFend', 'Human dORF Length', 'Human dORFseq', 'Homolog dORFstart', 'Homolog dORFend', 'Homolog dORF Length', 'Homolog dORFseq', 'Homolog dORF Length Difference', 'Homolog dORF Start Difference', 'Homolog dORF End Difference', 'dORF Water Length', 'dORF Water Score', 'dORF Water Identity', 'dORF Water Similarity', 'dORF Water Gaps', 'dORF vs 3UTR Water Identity Difference', 'dORF vs 3UTR Water Similarity Difference', 'dORF vs 3UTR Water Gaps Difference', 'Human CtrlStart', 'Human CtrlEnd', 'Human Ctrl Length', 'Human Ctrl Sequence', 'Ctrl Water Length', 'Ctrl Water Score', 'Ctrl Water Identity', 'Ctrl Water Similarity', 'Ctrl Water Gaps', 'Ctrl vs 3UTR Water Identity Difference', 'Ctrl vs 3UTR Water Similarity Difference', 'Ctrl vs 3UTR Water Gaps Difference']
ResultsDF_Ctrl = pd.DataFrame(columns=columns)

#Save Results DF as empty .csv with headings for columns
ResultsDF_Ctrl.to_csv('Human_and_Xtropicalis_ResultsDF_Ctrl_middle_water_included.csv')
print('Empty Results Dataframe .csv created')

#Create temporary subsets of results dataframes for each line in a loop
for i in range(len(ResultsDF)):
    #Create empty temp results dataframe
    tempResultsDF = pd.DataFrame(columns=columns)
    #Add row of results dataframe to temporary dataframe results
    tempResultsDF = pd.concat([tempResultsDF, ResultsDF.loc[[i]]])
    tempResultsDF = tempResultsDF.reset_index(drop=True)
    #Use human dORF length to gather location of control sequence from the middle of the human 3UTR of the same length as the dORF for that 3UTR
    HumandORFlength = tempResultsDF['Human dORF Length'][0]
    Human3UTRlength = tempResultsDF['Human 3UTR Length'][0]
    HumanCtrlstart = int(Human3UTRlength / 2)
    HumanCtrlend = HumanCtrlstart + HumandORFlength
    #Only continue loop if the end of the Ctrl seqs falls within the 3' UTR
    if HumanCtrlend < tempResultsDF['Human 3UTR Length'][0]:
        #Gather Ctrl sequences using locations as index for 3' UTR string and calculate the length and get start and end by adding 1 to the indexes - for homolog simply take the homolog 3UTR from the results dataframe
        HumanCtrlSeq = tempResultsDF['Human 3UTRseq'][0][HumanCtrlstart:HumanCtrlend]
        HumanCtrlstart = HumanCtrlstart + 1
        HumanCtrlend = HumanCtrlend
        HumanCtrlLen = len(HumanCtrlSeq)
        HomologCtrlSeq = tempResultsDF['Homolog 3UTRseq'][0]
        HomologCtrlLen = len(HomologCtrlSeq)
        #Write out Ctrl humand and homolog sequences into text files to be used with water tool in command line
        Human_Ctrl_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'XtroHuman_Ctrl_temp_seq.txt'), 'w')
        Human_Ctrl_temp_seq.write(HumanCtrlSeq)
        Human_Ctrl_temp_seq.close()
        Homolog_Ctrl_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'XtroHomolog_Ctrl_temp_seq.txt'), 'w')
        Homolog_Ctrl_temp_seq.write(HomologCtrlSeq)
        Homolog_Ctrl_temp_seq.close()

        #run water alignment with 2 temp 3UTR sequence files
        os.system("water XtroHuman_Ctrl_temp_seq.txt XtroHomolog_Ctrl_temp_seq.txt -gapopen 10.0 -gapextend 0.5 -outfile XtroCtrlTemp_Align.water >/dev/null 2>&1")

        #remove temporary sequence files - no longer needed - save space
        os.system("rm XtroHuman_Ctrl_temp_seq.txt")
        os.system("rm XtroHomolog_Ctrl_temp_seq.txt")

        #open water alignment file from path as string in python script
        with open('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis/XtroCtrlTemp_Align.water') as Ctrlwater_file:
            CtrlTemp_Align = Ctrlwater_file.read()
        
        #Take water alignment file string line by line and pull out lines matching the results wanted    
        for lines in CtrlTemp_Align.split("\n"):
            if "Length" in lines:
                CtrlTemp_Length = lines.strip()
        for lines in CtrlTemp_Align.split("\n"):
            if "Identity" in lines:
                CtrlTemp_Identity = lines.strip()
        for lines in CtrlTemp_Align.split("\n"):
            if "Similarity" in lines:
                CtrlTemp_Similarity = lines.strip()
        for lines in CtrlTemp_Align.split("\n"):
            if "Gaps" in lines:
                CtrlTemp_Gaps = lines.strip()
        for lines in CtrlTemp_Align.split("\n"):
            if "Score" in lines:
                CtrlTemp_Score = lines.strip()

        #use regular expressions to gather water alignment values from selected lines of alignment results file
        CtrlTemp_Length = re.search("[0-9]{1,}", CtrlTemp_Length)
        CtrlTemp_Length = CtrlTemp_Length.group(0)
        CtrlTemp_Length = int(CtrlTemp_Length)
        CtrlTemp_Score = re.search("[0-9]{1,}\.[0-9]", CtrlTemp_Score)
        CtrlTemp_Score = CtrlTemp_Score.group(0)
        CtrlTemp_Score = float(CtrlTemp_Score)
        CtrlTemp_Identity = re.search("[0-9]{1,}\.[0-9]", CtrlTemp_Identity)
        CtrlTemp_Identity = CtrlTemp_Identity.group(0)
        CtrlTemp_Identity = float(CtrlTemp_Identity)
        CtrlTemp_Similarity = re.search("[0-9]{1,}\.[0-9]", CtrlTemp_Similarity)
        CtrlTemp_Similarity = CtrlTemp_Similarity.group(0)
        CtrlTemp_Similarity = float(CtrlTemp_Similarity)
        CtrlTemp_Gaps = re.search("[0-9]{1,}\.[0-9]", CtrlTemp_Gaps)
        CtrlTemp_Gaps = CtrlTemp_Gaps.group(0)
        CtrlTemp_Gaps = float(CtrlTemp_Gaps)
        
        #remove Ctrl water alignment file from current directory using command line
        os.system("rm XtroCtrlTemp_Align.water")
        
        #if water alignment length is more or equal to the human or Ctrl length then add results to dataframe if not the whole of the human control sequence is not aligned
        if CtrlTemp_Length >= HumanCtrlLen:
            #add details about Ctrl sequence to temp results dataframe
            tempResultsDF['Human CtrlStart'][0] = HumanCtrlstart
            tempResultsDF['Human CtrlEnd'][0] = HumanCtrlend
            tempResultsDF['Human Ctrl Length'][0] = HumanCtrlLen
            tempResultsDF['Human Ctrl Sequence'][0] = HumanCtrlSeq
            
            #add Ctrl water alignment results to temp df - putting in the dORF loops means even if dORFs in same 3'UTRs as previous still add these details
            tempResultsDF['Ctrl Water Length'][0] = CtrlTemp_Length
            tempResultsDF['Ctrl Water Score'][0] = CtrlTemp_Score
            tempResultsDF['Ctrl Water Identity'][0] = CtrlTemp_Identity
            tempResultsDF['Ctrl Water Similarity'][0] = CtrlTemp_Similarity
            tempResultsDF['Ctrl Water Gaps'][0] = CtrlTemp_Gaps
            
            #Compare water alignments results for Ctrl compared to the 3' UTR between human and homolog
            threeUTRidentity = tempResultsDF['3UTR Water Identity'][0]
            tempIdentityDiff = CtrlTemp_Identity - threeUTRidentity
            tempIdentityDiff = format(tempIdentityDiff, '.2f')
            tempResultsDF['Ctrl vs 3UTR Water Identity Difference'] = [tempIdentityDiff]
            threeUTRsimilarity = tempResultsDF['3UTR Water Similarity'][0]
            tempSimilarityDiff = CtrlTemp_Similarity - threeUTRsimilarity
            tempSimilarityDiff = format(tempSimilarityDiff, '.2f')
            tempResultsDF['Ctrl vs 3UTR Water Similarity Difference'] = [tempSimilarityDiff]
            threeUTRgaps = tempResultsDF['3UTR Water Gaps'][0]
            tempGapsDiff = CtrlTemp_Gaps - threeUTRgaps
            tempGapsDiff = format(tempGapsDiff, '.2f')
            tempResultsDF['Ctrl vs 3UTR Water Gaps Difference'] = [tempGapsDiff]
            
            #write results of temporary results dataframe into created.csv results file
            tempResultsDF.to_csv('Human_and_Xtropicalis_ResultsDF_Ctrl_middle_water_included.csv', header=None, mode='a')
        else:
            continue
    else:
        continue
        
    #print rolling count of number of comparisons completed
    print('\r\033[K', end='')
    print("Number of results dataframe lines analysed", i, end='\r')
    
#print line to terminal to confirm analysis completed and dataframe saved
print('\n','Analysis Complete - Dataframe saved as csv')

