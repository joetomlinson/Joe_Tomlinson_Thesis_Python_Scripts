#Water analysis python script to take the 3’UTR and dORF output of human and a homolog [Change Homolog input script and output file for each homolog] to run comparisons between human and homolog in terms of smithwaterman alignments and also location and lengths of dORFs

#imported required modules
import pandas as pd
import os
import re
import warnings
from pandas.errors import SettingWithCopyWarning

#prevent warning showing up every time for the loops
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)


#load dataframes from .csv files and reset indexes
Human_rd = pd.read_csv("Human_RefSeq_dORF_3UTRs_remove_duplicate_3UTR.csv", index_col=0)
Human_rd = Human_rd.reset_index(drop=True)
Homologs_rd = pd.read_csv("Dmelanogaster_RefSeq_dORF_containing_3UTRs_remove_duplicate_3UTRs.csv", index_col=0)
Homologs_rd = Homologs_rd.reset_index(drop=True)

print('Human and Homolog Species dORF contatining 3UTRs imported')

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
    #Convert dORF IDs str into list of dORF sequence strings
    Human_rd['Human dORF ID'][i] = Human_rd['Human dORF ID'][i].strip('[')
    Human_rd['Human dORF ID'][i] = Human_rd['Human dORF ID'][i].strip(']')
    Human_rd['Human dORF ID'][i] = Human_rd['Human dORF ID'][i].strip("'")
    Human_rd['Human dORF ID'][i] = Human_rd['Human dORF ID'][i].split("', '")

#loop through whole homologs dataframe
for i in range(len(Homologs_rd)):
    #if only one dORF treated slightly differently but convert dORF location string to list of integers
    if Homologs_rd['Number of dORFs'][i] == 1:
        Homologs_rd['dORF Locations'][i] = Homologs_rd['dORF Locations'][i].strip('[')
        Homologs_rd['dORF Locations'][i] = Homologs_rd['dORF Locations'][i].strip(']')
        Homologs_rd['dORF Locations'][i] = Homologs_rd['dORF Locations'][i].split(', ')
        for n in range(len(Homologs_rd['dORF Locations'][i])):
            Homologs_rd['dORF Locations'][i][n] = int(Homologs_rd['dORF Locations'][i][n])
    #if more than one dORF convert dORF locations string to list of list of integers
    else:
       Homologs_rd['dORF Locations'][i] = Homologs_rd['dORF Locations'][i].split('], ')
       for x in range(len(Homologs_rd['dORF Locations'][i])):
           Homologs_rd['dORF Locations'][i][x] = Homologs_rd['dORF Locations'][i][x].strip('[')
           Homologs_rd['dORF Locations'][i][x] = Homologs_rd['dORF Locations'][i][x].strip(']')
           Homologs_rd['dORF Locations'][i][x] = Homologs_rd['dORF Locations'][i][x].split(', ')
           for n in range(len(Homologs_rd['dORF Locations'][i][x])):
               Homologs_rd['dORF Locations'][i][x][n] = int(Homologs_rd['dORF Locations'][i][x][n])
    #Convert dORF lengths string into list of int lengths of dORFs
    Homologs_rd['dORF Lengths'][i] = Homologs_rd['dORF Lengths'][i].strip('[')
    Homologs_rd['dORF Lengths'][i] = Homologs_rd['dORF Lengths'][i].strip(']')
    Homologs_rd['dORF Lengths'][i] = Homologs_rd['dORF Lengths'][i].split(", ")
    for y in range(len(Homologs_rd['dORF Lengths'][i])):
        Homologs_rd['dORF Lengths'][i][y] = int(Homologs_rd['dORF Lengths'][i][y])
    #Convert dORF sequences str into list of dORF sequence strings
    Homologs_rd['dORF Sequences'][i] = Homologs_rd['dORF Sequences'][i].strip('[')
    Homologs_rd['dORF Sequences'][i] = Homologs_rd['dORF Sequences'][i].strip(']')
    Homologs_rd['dORF Sequences'][i] = Homologs_rd['dORF Sequences'][i].strip("'")
    Homologs_rd['dORF Sequences'][i] = Homologs_rd['dORF Sequences'][i].split("', '")

print('Human and Homolog Dataframe reformatted')

#take list of human genes from dataframe and remove duplicates
HumanGenes = Homologs_rd['Human Gene Name'].tolist()
HumanGenes = list(set(HumanGenes))

print('Human Gene List Generated, number of human genes to try:', len(HumanGenes))

#Create results dataframe with list of columns as column headings
columns = ['Human Gene', 'Human Version', 'Human 3UTR Length', 'Human 3UTRseq', 'Homolog Species', 'Homolog Gene', 'Homolog Version', 'Homolog 3UTR Length', 'Homolog 3UTRseq', 'Homolog 3UTR Length Difference', '3UTR Water Length', '3UTR Water Score', '3UTR Water Identity', '3UTR Water Similarity', '3UTR Water Gaps', 'Human dORF ID', 'Human dORFstart', 'Human dORFend', 'Human dORF Length', 'Human dORFseq', 'Homolog dORFstart', 'Homolog dORFend', 'Homolog dORF Length', 'Homolog dORFseq', 'Homolog dORF Length Difference', 'Homolog dORF Start Difference', 'Homolog dORF End Difference', 'dORF Water Length', 'dORF Water Score', 'dORF Water Identity', 'dORF Water Similarity', 'dORF Water Gaps', 'dORF vs 3UTR Water Identity Difference', 'dORF vs 3UTR Water Similarity Difference', 'dORF vs 3UTR Water Gaps Difference']
ResultsDF = pd.DataFrame(columns=columns)

#Create temporary subsets of dataframes matching human gene
for i in range(len(HumanGenes)):
    tempHuman_rd = Human_rd.loc[Human_rd['Gene Name'] == HumanGenes[i]]
    tempHuman_rd = tempHuman_rd.reset_index(drop=True) #reset index do not keep from original dataframe
    tempHomologs_rd = Homologs_rd.loc[Homologs_rd['Human Gene Name'] == HumanGenes[i]]
    tempHomologs_rd = tempHomologs_rd.reset_index(drop=True)
    
#Create empty list to hold all homolog dORF seqs in homolog temp dataframe    
    tempHomologs_dORF_seqs = list()
    tempHomologs_dORF_ids = list()
#Loop through every dORF in homolog temp dataframe
    for d in range(len(tempHomologs_rd)):
        for z in range(0, tempHomologs_rd['Number of dORFs'][d]):
#Create variable to store each dORF sequence and add to list
            if tempHomologs_rd['Number of dORFs'][d] == 1:
                tempdORFseq = tempHomologs_rd['dORF Sequences'][d][0]
            else:
                tempdORFseq = tempHomologs_rd['dORF Sequences'][d][z]
            tempHomologs_dORF_seqs.append(tempdORFseq)
            dORF_ids = '>' + str(d) + '_' + str(z)
            tempHomologs_dORF_ids.append(dORF_ids)

#open text file in path
    Homologs_dORF_temp_seqs = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'Dmelanogaster_Homologs_dORF_temp_seqs.txt'), 'w')

#add each dORF in list to text file with ID in fasta format ready to run with water analysis
    for h in range(len(tempHomologs_dORF_seqs)):  
        Homologs_dORF_temp_seqs.write(tempHomologs_dORF_ids[h])
        Homologs_dORF_temp_seqs.write('\n')
        Homologs_dORF_temp_seqs.write(tempHomologs_dORF_seqs[h])
        Homologs_dORF_temp_seqs.write('\n')

#close the text file
    Homologs_dORF_temp_seqs.close()

#loop through each entry in temp human and homolog dataframes
    for a in range(len(tempHuman_rd)):
        for b in range(len(tempHomologs_rd)):
            
            #store the 3' UTR sequences from each df subset in loop as variables - ready for water alignment analysis
            Human_temp3UTR = tempHuman_rd['Three Prime Sequence'][a]
            Homologs_temp3UTR = tempHomologs_rd['Three Prime Sequence'][b]
            
            #Create temporary textfiles with 3UTR sequences from each dataframe subset in current directory - carried out in each loop for every comparison between 3UTRs
            Human_ThreeUTR_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'Dmelanogaster_Human_ThreeUTR_temp_seq.txt'), 'w')
            Human_ThreeUTR_temp_seq.write(Human_temp3UTR)
            Human_ThreeUTR_temp_seq.close()
            Homologs_ThreeUTR_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'Dmelanogaster_Homologs_ThreeUTR_temp_seq.txt'), 'w')
            Homologs_ThreeUTR_temp_seq.write(Homologs_temp3UTR)
            Homologs_ThreeUTR_temp_seq.close()

            #run water alignment with 2 temp 3UTR sequence files
            os.system("water Dmelanogaster_Human_ThreeUTR_temp_seq.txt Dmelanogaster_Homologs_ThreeUTR_temp_seq.txt -gapopen 10.0 -gapextend 0.5 -outfile Dmelanogaster_ThreeUTRTemp_Align.water >/dev/null 2>&1")

            #remove temporary sequence files - no longer needed - save space
            os.system("rm Dmelanogaster_Human_ThreeUTR_temp_seq.txt")
            os.system("rm Dmelanogaster_Homologs_ThreeUTR_temp_seq.txt")

            #open water alignment file from path as string in python script
            with open('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis/Dmelanogaster_ThreeUTRTemp_Align.water') as ThreeUTRwater_file:
                ThreeUTRTemp_Align = ThreeUTRwater_file.read()
            
            #Take water alignment file string line by line and pull out lines matching the results wanted    
            for lines in ThreeUTRTemp_Align.split("\n"):
                if "Length" in lines:
                    ThreeUTRTemp_Length = lines.strip()
            for lines in ThreeUTRTemp_Align.split("\n"):
                if "Identity" in lines:
                    ThreeUTRTemp_Identity = lines.strip()
            for lines in ThreeUTRTemp_Align.split("\n"):
                if "Similarity" in lines:
                    ThreeUTRTemp_Similarity = lines.strip()
            for lines in ThreeUTRTemp_Align.split("\n"):
                if "Gaps" in lines:
                    ThreeUTRTemp_Gaps = lines.strip()
            for lines in ThreeUTRTemp_Align.split("\n"):
                if "Score" in lines:
                    ThreeUTRTemp_Score = lines.strip()

            #use regular expressions to gather water alignment values from selected lines of alignment results file
            ThreeUTRTemp_Length = re.search("[0-9]{1,}", ThreeUTRTemp_Length)
            ThreeUTRTemp_Length = ThreeUTRTemp_Length.group(0)
            ThreeUTRTemp_Length = int(ThreeUTRTemp_Length)
            ThreeUTRTemp_Score = re.search("[0-9]{1,}\.[0-9]", ThreeUTRTemp_Score)
            ThreeUTRTemp_Score = ThreeUTRTemp_Score.group(0)
            ThreeUTRTemp_Score = float(ThreeUTRTemp_Score)
            ThreeUTRTemp_Identity = re.search("[0-9]{1,}\.[0-9]", ThreeUTRTemp_Identity)
            ThreeUTRTemp_Identity = ThreeUTRTemp_Identity.group(0)
            ThreeUTRTemp_Identity = float(ThreeUTRTemp_Identity)
            ThreeUTRTemp_Similarity = re.search("[0-9]{1,}\.[0-9]", ThreeUTRTemp_Similarity)
            ThreeUTRTemp_Similarity = ThreeUTRTemp_Similarity.group(0)
            ThreeUTRTemp_Similarity = float(ThreeUTRTemp_Similarity)
            ThreeUTRTemp_Gaps = re.search("[0-9]{1,}\.[0-9]", ThreeUTRTemp_Gaps)
            ThreeUTRTemp_Gaps = ThreeUTRTemp_Gaps.group(0)
            ThreeUTRTemp_Gaps = float(ThreeUTRTemp_Gaps)
            
            #remove 3UTR water alignment file from current directory using command line
            os.system("rm Dmelanogaster_ThreeUTRTemp_Align.water")

            #For the 3UTR seqs in the above loop from the subset df loop through each dORF in the 3UTR in human with each dORF in the 3' UTR in homologs for each loop - human 3UTR dORF looped with each homolog 3UTR in turn if matching human genes
            for x in range(0, tempHuman_rd['Number of dORFs'][a]):
                
                #Create new temporary dataframe to hold results of each loop
                tempResultsDF = pd.DataFrame(columns=columns)
                
                if tempHuman_rd['Number of dORFs'][a] == 1:
                    Human_tempdORFstart = tempHuman_rd['dORF Locations'][a][0]
                    Human_tempdORFend = tempHuman_rd['dORF Locations'][a][1]
                    Human_tempdORFlength = tempHuman_rd['dORF Lengths'][a][0]
                    Human_tempdORFseq = tempHuman_rd['dORF Sequences'][a][0]
                    tempHuman_dORF_ID = tempHuman_rd['Human dORF ID'][a][0]
                else:    
                    Human_tempdORFstart = tempHuman_rd['dORF Locations'][a][x][0]
                    Human_tempdORFend = tempHuman_rd['dORF Locations'][a][x][1]
                    Human_tempdORFlength = tempHuman_rd['dORF Lengths'][a][x]
                    Human_tempdORFseq = tempHuman_rd['dORF Sequences'][a][x]
                    tempHuman_dORF_ID = tempHuman_rd['Human dORF ID'][a][x]
                    
                #Create temporary textfiles with dORF sequences from each dataframe subset in current directory - carried out in each loop for every comparison between dORFs
                Human_dORF_temp_seq = open(os.path.join('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis', 'Dmelanogaster_Human_dORF_temp_seq.txt'), 'w')
                Human_dORF_temp_seq.write(Human_tempdORFseq)
                Human_dORF_temp_seq.close()

                #run water alignment with 2 temp dORF sequence files - homologs now have multiple dORF sequences
                os.system("water Dmelanogaster_Human_dORF_temp_seq.txt Dmelanogaster_Homologs_dORF_temp_seqs.txt -gapopen 10.0 -gapextend 0.5 -outfile Dmelanogaster_dORFTemp_Align.water >/dev/null 2>&1")

                #remove temporary sequence files - no longer needed - save space
                os.system("rm Dmelanogaster_Human_dORF_temp_seq.txt")

                #open water alignment file from path as string in python script
                with open('/users/payjt2/PhD/Homolog_3UTR_and_dORF_comparison_analysis/Emboss_Water_Analysis/Dmelanogaster_dORFTemp_Align.water') as dORFwater_file:
                    dORFTemp_Align = dORFwater_file.read()
                
                #Split dORF alignment file - will have multiple analysis results - split using #======== into list of strings
                split_dORFTemp_Align = dORFTemp_Align.split("#=======================================")
                #create empty list to store water results for each dORF human run against
                Seqs_dORFTemp_Align = list()
                #loop through list of water analysis results
                for Seqs in split_dORFTemp_Align:
                    #add section of results file to another list if it contains the phrase below
                    if "Aligned_sequences:" in Seqs:
                        Seqs_dORFTemp_Align.append(Seqs)

                #Create empty lists to store the water alignment results for each analysis done
                dORF_IDs = list()
                Lengths = list()
                Identities = list()
                Similarities = list()
                Gaps = list()
                Scores = list()

                #Loop through the listed water analysis results
                for c in range(len(Seqs_dORFTemp_Align)):
                #Take water alignment string line by line and pull out lines matching the results wanted
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "# 2: " in lines:
                            Temp_ID = lines.strip()
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "Length" in lines:
                            Temp_Length = lines.strip()
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "Identity" in lines:
                            Temp_Identity = lines.strip()
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "Similarity" in lines:
                            Temp_Similarity = lines.strip()
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "Gaps" in lines:
                            Temp_Gaps = lines.strip()
                    for lines in Seqs_dORFTemp_Align[c].split("\n"):
                        if "Score" in lines:
                            Temp_Score = lines.strip()

                #use regular expressions to gather water alignment values from selected lines of alignment results file - make sure appropriate variable eg int or float
                    Temp_ID = re.search("\:\ [0-9]{1,}\_[0-9]{1,}", Temp_ID)
                    Temp_ID = Temp_ID.group(0)
                    Temp_ID = Temp_ID.strip(':')
                    Temp_ID = Temp_ID.strip(' ')
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
                    
                    #Add the water alignment values collected to the empty lists made above
                    dORF_IDs.append(Temp_ID)
                    Lengths.append(Temp_Length)
                    Identities.append(Temp_Identity)
                    Similarities.append(Temp_Similarity)
                    Gaps.append(Temp_Gaps)
                    Scores.append(Temp_Score)
                
                #Calculate the largest similarity % in the list and store the index of this result - only take first appearance of largest - use index to pull out all results and dORF involved in analysis 
                max_similarity_index = Similarities.index(max(Similarities))
                
                homolog_dORF = dORF_IDs[max_similarity_index]
                homolog_dORF = homolog_dORF.split('_')
                j = int(homolog_dORF[0])
                k = int(homolog_dORF[1])
                
                #use command line to remove alignment results file - no longer needed - save space
                os.system("rm Dmelanogaster_dORFTemp_Align.water")
                
                if tempHomologs_rd['Number of dORFs'][j] == 1:
                    Homologs_tempdORFstart = tempHomologs_rd['dORF Locations'][j][0]
                    Homologs_tempdORFend = tempHomologs_rd['dORF Locations'][j][1]
                    Homologs_tempdORFlength = tempHomologs_rd['dORF Lengths'][j][0]
                    Homologs_tempdORFseq = tempHomologs_rd['dORF Sequences'][j][0]
                else:    
                    Homologs_tempdORFstart = tempHomologs_rd['dORF Locations'][j][k][0]
                    Homologs_tempdORFend = tempHomologs_rd['dORF Locations'][j][k][1]
                    Homologs_tempdORFlength = tempHomologs_rd['dORF Lengths'][j][k]
                    Homologs_tempdORFseq = tempHomologs_rd['dORF Sequences'][j][k]

                #add details about human and homolog 3UTRs for the loops analysed
                tempResultsDF['Human Gene'] = [tempHuman_rd['Gene Name'][a]]
                tempResultsDF['Human Version'] = [tempHuman_rd['Accession Version'][a]]
                tempResultsDF['Human 3UTR Length'] = [tempHuman_rd['Three Prime Length'][a]]
                tempResultsDF['Human 3UTRseq'] = [Human_temp3UTR]
                tempResultsDF['Homolog Species'] = [tempHomologs_rd['Species'][b]]
                tempResultsDF['Homolog Gene'] = [tempHomologs_rd['Gene Name'][b]]
                tempResultsDF['Homolog Version'] = [tempHomologs_rd['Accession Version'][b]]
                tempResultsDF['Homolog 3UTR Length'] = [tempHomologs_rd['Three Prime Length'][b]]
                tempResultsDF['Homolog 3UTRseq'] = [Homologs_temp3UTR]
                temp3UTRdiff = tempHomologs_rd['Three Prime Length'][b] - tempHuman_rd['Three Prime Length'][a]
                tempResultsDF['Homolog 3UTR Length Difference'] = [temp3UTRdiff]
                
                #add 3UTR water alignment results to temp df - putting in the dORF loops means even if dORFs in same 3'UTRs as previous still add these details
                tempResultsDF['3UTR Water Length'] = [ThreeUTRTemp_Length]
                tempResultsDF['3UTR Water Score'] = [ThreeUTRTemp_Score]
                tempResultsDF['3UTR Water Identity'] = [ThreeUTRTemp_Identity]
                tempResultsDF['3UTR Water Similarity'] = [ThreeUTRTemp_Similarity]
                tempResultsDF['3UTR Water Gaps'] = [ThreeUTRTemp_Gaps]
                
                #add details about human and homolog dORFs in the loops
                tempResultsDF['Human dORF ID'] = [tempHuman_dORF_ID]
                tempResultsDF['Human dORFstart'] = [Human_tempdORFstart]
                tempResultsDF['Human dORFend'] = [Human_tempdORFend]
                tempResultsDF['Human dORF Length'] = [Human_tempdORFlength]
                tempResultsDF['Human dORFseq'] = [Human_tempdORFseq]
                tempResultsDF['Homolog dORFstart'] = [Homologs_tempdORFstart]
                tempResultsDF['Homolog dORFend'] = [Homologs_tempdORFend]
                tempResultsDF['Homolog dORF Length'] = [Homologs_tempdORFlength]
                tempResultsDF['Homolog dORFseq'] = [Homologs_tempdORFseq]
                tempdORFlengthdiff = Homologs_tempdORFlength - Human_tempdORFlength
                tempResultsDF['Homolog dORF Length Difference'] = [tempdORFlengthdiff]
                tempdORFstartdiff = Homologs_tempdORFstart - Human_tempdORFstart
                tempResultsDF['Homolog dORF Start Difference'] = [tempdORFstartdiff]
                tempdORFenddiff = Homologs_tempdORFend - Human_tempdORFend
                tempResultsDF['Homolog dORF End Difference'] = [tempdORFenddiff]
                
                #add dORF water alignment results to temp df
                tempResultsDF['dORF Water Length'] = [Lengths[max_similarity_index]]
                tempResultsDF['dORF Water Score'] = [Scores[max_similarity_index]]
                tempResultsDF['dORF Water Identity'] = [Identities[max_similarity_index]]
                tempResultsDF['dORF Water Similarity'] = [Similarities[max_similarity_index]]
                tempResultsDF['dORF Water Gaps'] = [Gaps[max_similarity_index]]
                
                #Compare water alignments results for dORFs compared to the 3' UTR between human and homolog - look at difference between dORF and 3 UTR +ive suggests dORF has greater % when comparing
                tempIdentityDiff = Identities[max_similarity_index] - ThreeUTRTemp_Identity
                tempIdentityDiff = format(tempIdentityDiff, '.2f')
                tempResultsDF['dORF vs 3UTR Water Identity Difference'] = [tempIdentityDiff]
                tempSimilarityDiff = Similarities[max_similarity_index] - ThreeUTRTemp_Similarity
                tempSimilarityDiff = format(tempSimilarityDiff, '.2f')
                tempResultsDF['dORF vs 3UTR Water Similarity Difference'] = [tempSimilarityDiff]
                tempGapsDiff = Gaps[max_similarity_index] - ThreeUTRTemp_Gaps
                tempGapsDiff = format(tempGapsDiff, '.2f')
                tempResultsDF['dORF vs 3UTR Water Gaps Difference'] = [tempGapsDiff]
                
                #Add temporary dataframe results line to full results dataframe
                ResultsDF = pd.concat([ResultsDF, tempResultsDF])
                ResultsDF = ResultsDF.reset_index(drop=True)
                
                #print rolling count of number of comparisons completed
                print('\r\033[K', end='')
                print("Number of HumanGene entries run through water analysis:", i, end='\r')

#remove temporary sequence files - no longer needed - save space
    os.system("rm Dmelanogaster_Homologs_dORF_temp_seqs.txt")

#save dataframe as csv in current directory
ResultsDF.to_csv('Human_and_Dmelanogaster_ResultsDF.csv')
#print line to terminal to confirm analysis completed and dataframe saved
print('\n','Comparisons Complete - Dataframe saved as csv - Dataframe size: ', len(ResultsDF))

