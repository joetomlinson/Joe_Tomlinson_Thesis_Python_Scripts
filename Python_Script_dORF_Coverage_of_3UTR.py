# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:01:44 2023
@author: jtoml
"""

# Import pandas
import pandas as pd
  
# Import the excel file with 3 columns: transcript, ORF start, ORF stop from getorf output and call it xls_file
excel_file = pd.ExcelFile("C:\\Users\jtoml\Downloads\dORF_LengthsS3.xlsx") # input name S3 refers to shuffle repeat 3
  
# View the excel_file's sheet names
print(excel_file.sheet_names)
  
# Load the excel_file's Sheet1 as a dataframe
df = excel_file.parse('Sheet1')
print(df)

# Import Regular expressions
from re import search

# Take transcripts column from dataframe and convert to list, print that list then remove duplicates from the list and print a list without duplicates
transcripts_with_duplicates = list(df.transcript) # convert transcript column of df into a list
print(transcripts_with_duplicates) # print full list of transcript - including duplicates from those with multiple dORFs
transcripts = list(dict.fromkeys(transcripts_with_duplicates)) # remove the duplicated values within the list
print(transcripts) # print list of transcripts with the duplicates removed

##### Takes a list of transcripts and the original dataframe and changes to an alternate dataframe with each transcripts and the ORFstart and stop positions, 
##### if there are multiple dORFs per transcript the ORF start and stop values are stored in lists
List_for_df_update = list() # generate empty list for use in the loop to make new dataframe
for n in range(len(transcripts)): # generate loop to go through each of the transcripts from the transcripts list
    print(n) # shows how far through the transcripts the code is    
    ORF_T_Start = list() # generate empty list ready for dORF start values
    ORF_T_Stop = list() # generate empty list ready for dORF stop values
    for i in df.index: # generate loop to go through each row of df
        if search(transcripts[n], df['transcript'][i]): # if row of df matches the transcript the do the next steps
            temp_start = df['ORFstart'][i] # store matching row's dORF start value in temporary start
            temp_stop = df['ORFstop'][i] # store matching row's dORF stop value in temporary stop
            ORF_T_Start.append(temp_start) # add these values to the empty lists added previously allows multiple values for transcripts with multiple dORFs
            ORF_T_Stop.append(temp_stop)
    List_for_df_update.append((transcripts[n], ORF_T_Start, ORF_T_Stop)) # add these to the empty list from before the loop and add the transcript and multiple or single start and stop values

df_update = pd.DataFrame(data = List_for_df_update, columns=['transcript', 'ORFstart', 'ORFstop']) # use the previous list as data to populate a new dataset df_update with assigned columns and one row for each transcript
print(df_update) # print new dataframe

# Generate new dataset called df_range with two columns with 
# transcript and dORF ranges where if there are multiple dORFs the ranges are in a list
##### Slight change to previous removes the word range from the list of ranges - improved
List_for_df_range = list() # generate empty list for use in the loop to make new dataframe
for i in df_update.index: # generate loop to go through each row of df_update
    n = 0 # resets n to 0 for each loop through the rows of df_update
    dORF_range = list() # generate empty list ready for the start and stop of each dORF as a range
    while n < len(df_update['ORFstart'][i]): # generate loop to go through each dORF start and stop for each transcript
        temp_dORF_range = (df_update['ORFstart'][i][n], df_update['ORFstop'][i][n]) # temporary variable to store the dORF start and stop of each dORF
        dORF_range.append(temp_dORF_range) # add the dORF range from the temporary file to the empty or growing list
        n += 1 # continues looping through each dORF of a transcript
    List_for_df_range.append((df_update['transcript'][i], dORF_range)) # add dORF ranges with transcript ID as a list from the whole dataframe

df_range = pd.DataFrame(data = List_for_df_range, columns=['transcript', 'dORFrange']) # creates new dataframe with each transcript and the dORF ranges for each transcript as a list
print(df_range) # print new dataframe

#### Takes df_range dataframe and produce a new dataframe with each transcript with the
#### and the length of the overlapping dORFs included, also then adds a column for the total dORF length within the transcript
List_for_df_final = list() # generate empty list for use in the loop to make new dataframe
def overlap(x, y): # define overlap function
    return range(max(x[0], y[0]), min(x[-1], y[-1]) + 1) # overlap function returns the start and end of the overlapping ranges that overlap, multiple results if there are multiple overlaps in a list of ranges

for n in df_range.index: # loop to go through each row of the df_range dataframe
    dORF_length = list() # create empty list to store the lengths of dORFs both the single ones or overlapping ones
    ranges = df_range['dORFrange'][n] # ranges variable to store the list of dORF ranges for each transcript as it works through dataframe
    dORFstart, dORFstop = min([x[0] for x in ranges]), 0 # variables to store the start and stop of dORFs
    for i in ranges: # loop to go through the list of dORF ranges in ranges variable
        if i[0] == dORFstart: # if the dORF start is the same as the first range then dORF stop if next in list
            dORFstop = i[1]
    while dORFstart: # create loop to work through the dORF start and stops to identify overlaps and change dORF start and stop variables accordingly
        for _ in ranges:
            for i in ranges:
                if i[1] > dORFstop and overlap(i, [dORFstart, dORFstop]):
                    dORFstop = i[1]
        temp_dORF_length = dORFstop - dORFstart + 1 # calculates temporary length of single or overlapped dORFs for each overlap return
        dORF_length.append(temp_dORF_length) # add this temporary length value to the dORF length list
        try: # tests subsequent block of code for errors to try to prevent issues with code running
            dORFstart = min([x[0] for x in ranges if x[0] > dORFstop])
            for i in ranges:
                if i[0] == dORFstart:
                    dORFstop = i[1]
        except ValueError: # prevents error from causing larger issue instead stops this part of loop
            dORFstart = None
    List_for_df_final.append((df_update['transcript'][n], dORF_length)) # Add list of transcript with the list of dORF lengths ready to make new dataframe

df_final = pd.DataFrame(data = List_for_df_final, columns=['transcript', 'dORFlength']) # creates new dataframe with each transcript and the dORF lengths overlapped or not for each transcript as a list

df_final['Total_dORF_Length'] = "" # Create new empty column
n=0 # reset n as 0
while n < len(df_final): # loop through each row of the df_final dataframe to make new column with total dORF length
    if len(df_final['dORFlength'][n]) == 1: # if there is only one dORF length value the value in the total column is the same
        df_final['Total_dORF_Length'][n] = df_final['dORFlength'][n]
    else: # if there is more than one value in dORF length then value in new total column is the sum of the dORF length column for each transcript
        df_final['Total_dORF_Length'][n] = sum(df_final['dORFlength'][n])
    n += 1 # continue loop move to the next row
print(df_final) # print the dataframe with the new column

df_final.to_excel(r"C:\\Users\jtoml\Downloads\dORF_Lengths_python_outS3.xlsx", index=False) # export the dataframe as excel file, in output name S3 refers to shuffle repeat 3