# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:55:38 2024

@author: jtoml
"""

#Python Script to take list of SRA datasets and trim adapters and based on quality with fastqc to check results, then align these processed datasets against the human genome and filter out the number of RP reads in each dataset aligned to dORFs, 3' UTRs and genes of interest.

#import required modules
import os
import re
import pandas as pd
import warnings
from pandas.errors import SettingWithCopyWarning

#prevent warning showing up every time for the loops
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

#move to desired directory to carry out analysis
os.system("cd ~/PhD/RP_dataset_dORF_alignment_analysis")

#make sure linux modules are loaded
os.system("module load sra-tools")

print('Moved to directory for analysis and loaded sra-tools')

#list of RP datasets with SRA accession to run through script
SRA_Datasets = list(['SRR10850865', 'SRR10850866', 'SRR10850867', 'SRR10850868', 'SRR10850869', 'SRR10850870', 'SRR18592261', 'SRR18592262', 'SRR18592263', 'SRR18592264', 'SRR18592265', 'SRR18592266', 'SRR18592267', 'SRR18592268', 'SRR18592269', 'SRR18592270', 'SRR18592271', 'SRR18592272', 'SRR18592273', 'SRR18592274', 'SRR18592275', 'SRR18592276', 'SRR18592277', 'SRR18592278', 'SRR18592279'])

#Generate empty results table to gather results with column titles in list
columns = ['Dataset', 'ACTB RP Alignments to Gene', 'ACTB RP Alignments to 3UTR', 'TUBB RP Alignments to Gene', 'TUBB RP Alignments to 3UTR', 'GAPDH RP Alignments to Gene', 'GAPDH RP Alignments to 3UTR', 'dORF_10_NM_005751.5 RP Alignments to Gene', 'dORF_10_NM_005751.5 RP Alignments to 3UTR', 'dORF_10_NM_005751.5 RP Alignments to dORF', 'dORF_16_NM_001137675.4 RP Alignments to Gene', 'dORF_16_NM_001137675.4 RP Alignments to 3UTR', 'dORF_16_NM_001137675.4 RP Alignments to dORF', 'dORF_27_NM_001177387.1 RP Alignments to Gene', 'dORF_27_NM_001177387.1 RP Alignments to 3UTR', 'dORF_27_NM_001177387.1 RP Alignments to dORF', 'dORF_39_NM_001177387.1 RP Alignments to Gene', 'dORF_39_NM_001177387.1 RP Alignments to 3UTR', 'dORF_39_NM_001177387.1 RP Alignments to dORF', 'dORF_21_NM_001382309.1 RP Alignments to Gene', 'dORF_21_NM_001382309.1 RP Alignments to 3UTR', 'dORF_21_NM_001382309.1 RP Alignments to dORF', 'dORF_14_NM_006317.5 RP Alignments to Gene', 'dORF_14_NM_006317.5 RP Alignments to 3UTR', 'dORF_14_NM_006317.5 RP Alignments to dORF', 'dORF_17_NM_006317.5 RP Alignments to Gene', 'dORF_17_NM_006317.5 RP Alignments to 3UTR', 'dORF_17_NM_006317.5 RP Alignments to dORF', 'dORF_4_NM_001370402.1 RP Alignments to Gene', 'dORF_4_NM_001370402.1 RP Alignments to 3UTR', 'dORF_4_NM_001370402.1 RP Alignments to dORF', 'dORF_18_NM_001300905.2 RP Alignments to Gene', 'dORF_18_NM_001300905.2 RP Alignments to 3UTR', 'dORF_18_NM_001300905.2 RP Alignments to dORF', 'dORF_3_NM_001300905.2 RP Alignments to Gene', 'dORF_3_NM_001300905.2 RP Alignments to 3UTR', 'dORF_3_NM_001300905.2 RP Alignments to dORF', 'dORF_4_NM_001300905.2 RP Alignments to Gene', 'dORF_4_NM_001300905.2 RP Alignments to 3UTR', 'dORF_4_NM_001300905.2 RP Alignments to dORF', 'dORF_6_XM_054341607.1 RP Alignments to Gene', 'dORF_6_XM_054341607.1 RP Alignments to 3UTR', 'dORF_6_XM_054341607.1 RP Alignments to dORF', 'dORF_64_XM_024453667.2 RP Alignments to Gene', 'dORF_64_XM_024453667.2 RP Alignments to 3UTR', 'dORF_64_XM_024453667.2 RP Alignments to dORF', 'dORF_65_XM_024453667.2 RP Alignments to Gene', 'dORF_65_XM_024453667.2 RP Alignments to 3UTR', 'dORF_65_XM_024453667.2 RP Alignments to dORF', 'dORF_8_NM_001200.4 RP Alignments to Gene', 'dORF_8_NM_001200.4 RP Alignments to 3UTR', 'dORF_8_NM_001200.4 RP Alignments to dORF', 'dORF_6_NM_001204.7 RP Alignments to Gene', 'dORF_6_NM_001204.7 RP Alignments to 3UTR', 'dORF_6_NM_001204.7 RP Alignments to dORF', 'dORF_16_XM_047449149.1 RP Alignments to Gene', 'dORF_16_XM_047449149.1 RP Alignments to 3UTR', 'dORF_16_XM_047449149.1 RP Alignments to dORF', 'dORF_56_XM_054328991.1 RP Alignments to Gene', 'dORF_56_XM_054328991.1 RP Alignments to 3UTR', 'dORF_56_XM_054328991.1 RP Alignments to dORF', 'dORF_2_NM_001321690.1 RP Alignments to Gene', 'dORF_2_NM_001321690.1 RP Alignments to 3UTR', 'dORF_2_NM_001321690.1 RP Alignments to dORF', 'dORF_16_NM_000719.7 RP Alignments to Gene', 'dORF_16_NM_000719.7 RP Alignments to 3UTR', 'dORF_16_NM_000719.7 RP Alignments to dORF', 'dORF_4_XM_005252591.4 RP Alignments to Gene', 'dORF_4_XM_005252591.4 RP Alignments to 3UTR', 'dORF_4_XM_005252591.4 RP Alignments to dORF', 'dORF_1_NM_172128.3 RP Alignments to Gene', 'dORF_1_NM_172128.3 RP Alignments to 3UTR', 'dORF_1_NM_172128.3 RP Alignments to dORF', 'dORF_5_XM_047426959.1 RP Alignments to Gene', 'dORF_5_XM_047426959.1 RP Alignments to 3UTR', 'dORF_5_XM_047426959.1 RP Alignments to dORF', 'dORF_2_XM_047426959.1 RP Alignments to Gene', 'dORF_2_XM_047426959.1 RP Alignments to 3UTR', 'dORF_2_XM_047426959.1 RP Alignments to dORF', 'dORF_70_XM_054379598.1 RP Alignments to Gene', 'dORF_70_XM_054379598.1 RP Alignments to 3UTR', 'dORF_70_XM_054379598.1 RP Alignments to dORF', 'dORF_7_NM_001300829.2 RP Alignments to Gene', 'dORF_7_NM_001300829.2 RP Alignments to 3UTR', 'dORF_7_NM_001300829.2 RP Alignments to dORF', 'dORF_10_XM_047447764.1 RP Alignments to Gene', 'dORF_10_XM_047447764.1 RP Alignments to 3UTR', 'dORF_10_XM_047447764.1 RP Alignments to dORF', 'dORF_7_NM_001204062.2 RP Alignments to Gene', 'dORF_7_NM_001204062.2 RP Alignments to 3UTR', 'dORF_7_NM_001204062.2 RP Alignments to dORF', 'dORF_11_XM_017018086.2 RP Alignments to Gene', 'dORF_11_XM_017018086.2 RP Alignments to 3UTR', 'dORF_11_XM_017018086.2 RP Alignments to dORF', 'dORF_12_XM_017018086.2 RP Alignments to Gene', 'dORF_12_XM_017018086.2 RP Alignments to 3UTR', 'dORF_12_XM_017018086.2 RP Alignments to dORF', 'dORF_16_XM_054359700.1 RP Alignments to Gene', 'dORF_16_XM_054359700.1 RP Alignments to 3UTR', 'dORF_16_XM_054359700.1 RP Alignments to dORF', 'dORF_6_NM_001288715.1 RP Alignments to Gene', 'dORF_6_NM_001288715.1 RP Alignments to 3UTR', 'dORF_6_NM_001288715.1 RP Alignments to dORF', 'dORF_17_NM_080759.6 RP Alignments to Gene', 'dORF_17_NM_080759.6 RP Alignments to 3UTR', 'dORF_17_NM_080759.6 RP Alignments to dORF', 'dORF_16_XM_054313463.1 RP Alignments to Gene', 'dORF_16_XM_054313463.1 RP Alignments to 3UTR', 'dORF_16_XM_054313463.1 RP Alignments to dORF', 'dORF_32_XM_024453007.2 RP Alignments to Gene', 'dORF_32_XM_024453007.2 RP Alignments to 3UTR', 'dORF_32_XM_024453007.2 RP Alignments to dORF', 'dORF_18_NM_001364157.2 RP Alignments to Gene', 'dORF_18_NM_001364157.2 RP Alignments to 3UTR', 'dORF_18_NM_001364157.2 RP Alignments to dORF', 'dORF_7_NM_001964.3 RP Alignments to Gene', 'dORF_7_NM_001964.3 RP Alignments to 3UTR', 'dORF_7_NM_001964.3 RP Alignments to dORF', 'dORF_45_NM_001243513.1 RP Alignments to Gene', 'dORF_45_NM_001243513.1 RP Alignments to 3UTR', 'dORF_45_NM_001243513.1 RP Alignments to dORF', 'dORF_22_NM_016605.3 RP Alignments to Gene', 'dORF_22_NM_016605.3 RP Alignments to 3UTR', 'dORF_22_NM_016605.3 RP Alignments to dORF', 'dORF_23_NM_203301.4 RP Alignments to Gene', 'dORF_23_NM_203301.4 RP Alignments to 3UTR', 'dORF_23_NM_203301.4 RP Alignments to dORF', 'dORF_14_NM_003868.3 RP Alignments to Gene', 'dORF_14_NM_003868.3 RP Alignments to 3UTR', 'dORF_14_NM_003868.3 RP Alignments to dORF', 'dORF_30_NM_001346114.2 RP Alignments to Gene', 'dORF_30_NM_001346114.2 RP Alignments to 3UTR', 'dORF_30_NM_001346114.2 RP Alignments to dORF', 'dORF_3_XM_054339463.1 RP Alignments to Gene', 'dORF_3_XM_054339463.1 RP Alignments to 3UTR', 'dORF_3_XM_054339463.1 RP Alignments to dORF', 'dORF_6_NM_138426.4 RP Alignments to Gene', 'dORF_6_NM_138426.4 RP Alignments to 3UTR', 'dORF_6_NM_138426.4 RP Alignments to dORF', 'dORF_9_NM_022130.4 RP Alignments to Gene', 'dORF_9_NM_022130.4 RP Alignments to 3UTR', 'dORF_9_NM_022130.4 RP Alignments to dORF', 'dORF_2_NM_007325.5 RP Alignments to Gene', 'dORF_2_NM_007325.5 RP Alignments to 3UTR', 'dORF_2_NM_007325.5 RP Alignments to dORF', 'dORF_2_NM_024503.5 RP Alignments to Gene', 'dORF_2_NM_024503.5 RP Alignments to 3UTR', 'dORF_2_NM_024503.5 RP Alignments to dORF', 'dORF_11_XM_054358064.1 RP Alignments to Gene', 'dORF_11_XM_054358064.1 RP Alignments to 3UTR', 'dORF_11_XM_054358064.1 RP Alignments to dORF', 'dORF_17_NM_001363050.1 RP Alignments to Gene', 'dORF_17_NM_001363050.1 RP Alignments to 3UTR', 'dORF_17_NM_001363050.1 RP Alignments to dORF', 'dORF_20_NM_001363050.1 RP Alignments to Gene', 'dORF_20_NM_001363050.1 RP Alignments to 3UTR', 'dORF_20_NM_001363050.1 RP Alignments to dORF', 'dORF_19_NM_002253.4 RP Alignments to Gene', 'dORF_19_NM_002253.4 RP Alignments to 3UTR', 'dORF_19_NM_002253.4 RP Alignments to dORF', 'dORF_13_XM_011532501.3 RP Alignments to Gene', 'dORF_13_XM_011532501.3 RP Alignments to 3UTR', 'dORF_13_XM_011532501.3 RP Alignments to dORF', 'dORF_2_XM_047430302.1 RP Alignments to Gene', 'dORF_2_XM_047430302.1 RP Alignments to 3UTR', 'dORF_2_XM_047430302.1 RP Alignments to dORF', 'dORF_54_XM_054331699.1 RP Alignments to Gene', 'dORF_54_XM_054331699.1 RP Alignments to 3UTR', 'dORF_54_XM_054331699.1 RP Alignments to dORF', 'dORF_11_NM_004631.5 RP Alignments to Gene', 'dORF_11_NM_004631.5 RP Alignments to 3UTR', 'dORF_11_NM_004631.5 RP Alignments to dORF', 'dORF_4_NM_004631.5 RP Alignments to Gene', 'dORF_4_NM_004631.5 RP Alignments to 3UTR', 'dORF_4_NM_004631.5 RP Alignments to dORF', 'dORF_3_NM_001350216.3 RP Alignments to Gene', 'dORF_3_NM_001350216.3 RP Alignments to 3UTR', 'dORF_3_NM_001350216.3 RP Alignments to dORF', 'dORF_5_NM_015578.4 RP Alignments to Gene', 'dORF_5_NM_015578.4 RP Alignments to 3UTR', 'dORF_5_NM_015578.4 RP Alignments to dORF', 'dORF_8_NM_001351625.3 RP Alignments to Gene', 'dORF_8_NM_001351625.3 RP Alignments to 3UTR', 'dORF_8_NM_001351625.3 RP Alignments to dORF', 'dORF_8_XM_054346559.1 RP Alignments to Gene', 'dORF_8_XM_054346559.1 RP Alignments to 3UTR', 'dORF_8_XM_054346559.1 RP Alignments to dORF', 'dORF_8_NM_007358.4 RP Alignments to Gene', 'dORF_8_NM_007358.4 RP Alignments to 3UTR', 'dORF_8_NM_007358.4 RP Alignments to dORF', 'dORF_9_NM_007358.4 RP Alignments to Gene', 'dORF_9_NM_007358.4 RP Alignments to 3UTR', 'dORF_9_NM_007358.4 RP Alignments to dORF', 'dORF_45_NM_001329851.3 RP Alignments to Gene', 'dORF_45_NM_001329851.3 RP Alignments to 3UTR', 'dORF_45_NM_001329851.3 RP Alignments to dORF', 'dORF_7_NM_021076.4 RP Alignments to Gene', 'dORF_7_NM_021076.4 RP Alignments to 3UTR', 'dORF_7_NM_021076.4 RP Alignments to dORF', 'dORF_7_NM_001135659.3 RP Alignments to Gene', 'dORF_7_NM_001135659.3 RP Alignments to 3UTR', 'dORF_7_NM_001135659.3 RP Alignments to dORF', 'dORF_2_XM_047418839.1 RP Alignments to Gene', 'dORF_2_XM_047418839.1 RP Alignments to 3UTR', 'dORF_2_XM_047418839.1 RP Alignments to dORF', 'dORF_12_XM_047430276.1 RP Alignments to Gene', 'dORF_12_XM_047430276.1 RP Alignments to 3UTR', 'dORF_12_XM_047430276.1 RP Alignments to dORF', 'dORF_14_XM_047430276.1 RP Alignments to Gene', 'dORF_14_XM_047430276.1 RP Alignments to 3UTR', 'dORF_14_XM_047430276.1 RP Alignments to dORF', 'dORF_4_XM_011535918.4 RP Alignments to Gene', 'dORF_4_XM_011535918.4 RP Alignments to 3UTR', 'dORF_4_XM_011535918.4 RP Alignments to dORF', 'dORF_22_NM_001382323.2 RP Alignments to Gene', 'dORF_22_NM_001382323.2 RP Alignments to 3UTR', 'dORF_22_NM_001382323.2 RP Alignments to dORF', 'dORF_23_NM_001382323.2 RP Alignments to Gene', 'dORF_23_NM_001382323.2 RP Alignments to 3UTR', 'dORF_23_NM_001382323.2 RP Alignments to dORF', 'dORF_1_XM_054331366.1 RP Alignments to Gene', 'dORF_1_XM_054331366.1 RP Alignments to 3UTR', 'dORF_1_XM_054331366.1 RP Alignments to dORF', 'dORF_2_XM_017001874.2 RP Alignments to Gene', 'dORF_2_XM_017001874.2 RP Alignments to 3UTR', 'dORF_2_XM_017001874.2 RP Alignments to dORF', 'dORF_32_NM_130393.3 RP Alignments to Gene', 'dORF_32_NM_130393.3 RP Alignments to 3UTR', 'dORF_32_NM_130393.3 RP Alignments to dORF', 'dORF_19_NM_001020658.2 RP Alignments to Gene', 'dORF_19_NM_001020658.2 RP Alignments to 3UTR', 'dORF_19_NM_001020658.2 RP Alignments to dORF', 'dORF_31_NM_206854.3 RP Alignments to Gene', 'dORF_31_NM_206854.3 RP Alignments to 3UTR', 'dORF_31_NM_206854.3 RP Alignments to dORF', 'dORF_5_NM_001378107.1 RP Alignments to Gene', 'dORF_5_NM_001378107.1 RP Alignments to 3UTR', 'dORF_5_NM_001378107.1 RP Alignments to dORF', 'dORF_26_NM_001100588.3 RP Alignments to Gene', 'dORF_26_NM_001100588.3 RP Alignments to 3UTR', 'dORF_26_NM_001100588.3 RP Alignments to dORF', 'dORF_22_XM_054361422.1 RP Alignments to Gene', 'dORF_22_XM_054361422.1 RP Alignments to 3UTR', 'dORF_22_XM_054361422.1 RP Alignments to dORF', 'dORF_23_XM_054361422.1 RP Alignments to Gene', 'dORF_23_XM_054361422.1 RP Alignments to 3UTR', 'dORF_23_XM_054361422.1 RP Alignments to dORF', 'dORF_19_NM_175634.3 RP Alignments to Gene', 'dORF_19_NM_175634.3 RP Alignments to 3UTR', 'dORF_19_NM_175634.3 RP Alignments to dORF', 'dORF_23_NM_175634.3 RP Alignments to Gene', 'dORF_23_NM_175634.3 RP Alignments to 3UTR', 'dORF_23_NM_175634.3 RP Alignments to dORF', 'dORF_3_NM_175634.3 RP Alignments to Gene', 'dORF_3_NM_175634.3 RP Alignments to 3UTR', 'dORF_3_NM_175634.3 RP Alignments to dORF', 'dORF_5_NM_175634.3 RP Alignments to Gene', 'dORF_5_NM_175634.3 RP Alignments to 3UTR', 'dORF_5_NM_175634.3 RP Alignments to dORF', 'dORF_1_XM_054375659.1 RP Alignments to Gene', 'dORF_1_XM_054375659.1 RP Alignments to 3UTR', 'dORF_1_XM_054375659.1 RP Alignments to dORF', 'dORF_9_NM_002971.6 RP Alignments to Gene', 'dORF_9_NM_002971.6 RP Alignments to 3UTR', 'dORF_9_NM_002971.6 RP Alignments to dORF', 'dORF_20_NM_001172509.2 RP Alignments to Gene', 'dORF_20_NM_001172509.2 RP Alignments to 3UTR', 'dORF_20_NM_001172509.2 RP Alignments to dORF', 'dORF_29_NM_001358351.3 RP Alignments to Gene', 'dORF_29_NM_001358351.3 RP Alignments to 3UTR', 'dORF_29_NM_001358351.3 RP Alignments to dORF', 'dORF_11_NM_015559.3 RP Alignments to Gene', 'dORF_11_NM_015559.3 RP Alignments to 3UTR', 'dORF_11_NM_015559.3 RP Alignments to dORF', 'dORF_7_NM_001372044.2 RP Alignments to Gene', 'dORF_7_NM_001372044.2 RP Alignments to 3UTR', 'dORF_7_NM_001372044.2 RP Alignments to dORF', 'dORF_5_XM_047426145.1 RP Alignments to Gene', 'dORF_5_XM_047426145.1 RP Alignments to 3UTR', 'dORF_5_XM_047426145.1 RP Alignments to dORF', 'dORF_6_XM_047426145.1 RP Alignments to Gene', 'dORF_6_XM_047426145.1 RP Alignments to 3UTR', 'dORF_6_XM_047426145.1 RP Alignments to dORF', 'dORF_1_NM_005413.4 RP Alignments to Gene', 'dORF_1_NM_005413.4 RP Alignments to 3UTR', 'dORF_1_NM_005413.4 RP Alignments to dORF', 'dORF_21_NM_001304421.2 RP Alignments to Gene', 'dORF_21_NM_001304421.2 RP Alignments to 3UTR', 'dORF_21_NM_001304421.2 RP Alignments to dORF', 'dORF_2_XM_054321826.1 RP Alignments to Gene', 'dORF_2_XM_054321826.1 RP Alignments to 3UTR', 'dORF_2_XM_054321826.1 RP Alignments to dORF', 'dORF_1_XM_054365273.1 RP Alignments to Gene', 'dORF_1_XM_054365273.1 RP Alignments to 3UTR', 'dORF_1_XM_054365273.1 RP Alignments to dORF', 'dORF_5_NM_001172712.1 RP Alignments to Gene', 'dORF_5_NM_001172712.1 RP Alignments to 3UTR', 'dORF_5_NM_001172712.1 RP Alignments to dORF', 'dORF_6_NM_001309444.2 RP Alignments to Gene', 'dORF_6_NM_001309444.2 RP Alignments to 3UTR', 'dORF_6_NM_001309444.2 RP Alignments to dORF', 'dORF_4_XM_047436525.1 RP Alignments to Gene', 'dORF_4_XM_047436525.1 RP Alignments to 3UTR', 'dORF_4_XM_047436525.1 RP Alignments to dORF', 'dORF_11_NM_001169117.2 RP Alignments to Gene', 'dORF_11_NM_001169117.2 RP Alignments to 3UTR', 'dORF_11_NM_001169117.2 RP Alignments to dORF', 'dORF_4_NM_182910.2 RP Alignments to Gene', 'dORF_4_NM_182910.2 RP Alignments to 3UTR', 'dORF_4_NM_182910.2 RP Alignments to dORF', 'dORF_1_XM_054356666.1 RP Alignments to Gene', 'dORF_1_XM_054356666.1 RP Alignments to 3UTR', 'dORF_1_XM_054356666.1 RP Alignments to dORF', 'dORF_63_XM_017024430.3 RP Alignments to Gene', 'dORF_63_XM_017024430.3 RP Alignments to 3UTR', 'dORF_63_XM_017024430.3 RP Alignments to dORF', 'dORF_1_XM_054325883.1 RP Alignments to Gene', 'dORF_1_XM_054325883.1 RP Alignments to 3UTR', 'dORF_1_XM_054325883.1 RP Alignments to dORF', 'dORF_2_XM_054325883.1 RP Alignments to Gene', 'dORF_2_XM_054325883.1 RP Alignments to 3UTR', 'dORF_2_XM_054325883.1 RP Alignments to dORF', 'dORF_1_XM_054325887.1 RP Alignments to Gene', 'dORF_1_XM_054325887.1 RP Alignments to 3UTR', 'dORF_1_XM_054325887.1 RP Alignments to dORF', 'dORF_26_NM_001098816.3 RP Alignments to Gene', 'dORF_26_NM_001098816.3 RP Alignments to 3UTR', 'dORF_26_NM_001098816.3 RP Alignments to dORF', 'dORF_6_XM_047436271.1 RP Alignments to Gene', 'dORF_6_XM_047436271.1 RP Alignments to 3UTR', 'dORF_6_XM_047436271.1 RP Alignments to dORF', 'dORF_13_XM_011518972.4 RP Alignments to Gene', 'dORF_13_XM_011518972.4 RP Alignments to 3UTR', 'dORF_13_XM_011518972.4 RP Alignments to dORF', 'dORF_20_NM_001195071.1 RP Alignments to Gene', 'dORF_20_NM_001195071.1 RP Alignments to 3UTR', 'dORF_20_NM_001195071.1 RP Alignments to dORF', 'dORF_17_XM_054316816.1 RP Alignments to Gene', 'dORF_17_XM_054316816.1 RP Alignments to 3UTR', 'dORF_17_XM_054316816.1 RP Alignments to dORF', 'dORF_14_NM_001080430.4 RP Alignments to Gene', 'dORF_14_NM_001080430.4 RP Alignments to 3UTR', 'dORF_14_NM_001080430.4 RP Alignments to dORF', 'dORF_23_NM_014112.5 RP Alignments to Gene', 'dORF_23_NM_014112.5 RP Alignments to 3UTR', 'dORF_23_NM_014112.5 RP Alignments to dORF', 'dORF_21_XM_054347743.1 RP Alignments to Gene', 'dORF_21_XM_054347743.1 RP Alignments to 3UTR', 'dORF_21_XM_054347743.1 RP Alignments to dORF', 'dORF_9_XM_047416346.1 RP Alignments to Gene', 'dORF_9_XM_047416346.1 RP Alignments to 3UTR', 'dORF_9_XM_047416346.1 RP Alignments to dORF', 'dORF_3_XM_047419923.1 RP Alignments to Gene', 'dORF_3_XM_047419923.1 RP Alignments to 3UTR', 'dORF_3_XM_047419923.1 RP Alignments to dORF', 'dORF_19_NM_017590.6 RP Alignments to Gene', 'dORF_19_NM_017590.6 RP Alignments to 3UTR', 'dORF_19_NM_017590.6 RP Alignments to dORF', 'dORF_17_NM_001278244.1 RP Alignments to Gene', 'dORF_17_NM_001278244.1 RP Alignments to 3UTR', 'dORF_17_NM_001278244.1 RP Alignments to dORF', 'dORF_38_XM_047434169.1 RP Alignments to Gene', 'dORF_38_XM_047434169.1 RP Alignments to 3UTR', 'dORF_38_XM_047434169.1 RP Alignments to dORF', 'dORF_28_XM_017027246.3 RP Alignments to Gene', 'dORF_28_XM_017027246.3 RP Alignments to 3UTR', 'dORF_28_XM_017027246.3 RP Alignments to dORF', 'dORF_23_XM_011526270.4 RP Alignments to Gene', 'dORF_23_XM_011526270.4 RP Alignments to 3UTR', 'dORF_23_XM_011526270.4 RP Alignments to dORF', 'dORF_75_XM_011521543.4 RP Alignments to Gene', 'dORF_75_XM_011521543.4 RP Alignments to 3UTR', 'dORF_75_XM_011521543.4 RP Alignments to dORF']
ResultsDF = pd.DataFrame(columns=columns)

#Save Results DF as empty .csv with headings for columns
ResultsDF.to_csv('Part6_RP_alignment_gene_3UTR_HAS_dORF_results.csv') ###modify name
print('Empty Results Dataframe .csv created')

#list of regions to look at RP read alignment
Housekeeping_Region_list = list([['ACTB', '-', '7:5527168-5530581', '7:5527167-5527727'], ['TUBB', '+', '6:30720372-30725402', '6:30724417-30725402'], ['GAPDH', '+', '12:6534537-6538351', '12:6538190-6538351']])
dORF_Region_List_HAS = list([['dORF_10_NM_005751.5', '+', '7:91940882-92110623', '7:92110179-92110653', '7:92110618-92110623'], ['dORF_16_NM_001137675.4', '+', '16:71845996-71857308', '16:71851830-71857308', '16:71853642-71853683'], ['dORF_27_NM_001177387.1', '+', '3:63863164-64003442', '3:63999579-64003442', '3:64002117-64002200'], ['dORF_39_NM_001177387.1', '+', '3:63863164-64003442', '3:63999579-64003442', '3:64003362-64003412'], ['dORF_21_NM_001382309.1', '-', '17:44191825-44199864', '17:44191824-44194242', '17:44191940-44191954'], ['dORF_14_NM_006317.5', '+', '5:17216843-17276814', '5:17275920-17276814', '5:17276650-17276808'], ['dORF_17_NM_006317.5', '+', '5:17216843-17276814', '5:17275920-17276814', '5:17276777-17276782'], ['dORF_4_NM_001370402.1', '-', '7:73440426-73522273', '7:73440425-73442175', '7:73441400-73441420'], ['dORF_18_NM_001300905.2', '-', '12:56595616-56638298', '12:56595615-56598597', '12:56596123-56596296'], ['dORF_3_NM_001300905.2', '-', '12:56595616-56638298', '12:56595615-56598597', '12:56598181-56598306'], ['dORF_4_NM_001300905.2', '-', '12:56595616-56638298', '12:56595615-56598597', '12:56598191-56598199'], ['dORF_6_XM_054341607.1', '-', '2:159315332-159712422', '2:159318998-159320244', '2:159319499-159319564'], ['dORF_64_XM_024453667.2', '+', '3:107522982-107811319', '3:107805547-107811319', '3:107811014-107811151'], ['dORF_65_XM_024453667.2', '+', '3:107522982-107811319', '3:107805547-107811319', '3:107811047-107811151'], ['dORF_8_NM_001200.4', '+', '20:6767706-6780226', '20:6779109-6780226', '20:6779396-6779404'], ['dORF_6_NM_001204.7', '+', '2:202376347-202567729', '2:202559966-202567729', '2:202560227-202560256'], ['dORF_16_XM_047449149.1', '+', '3:49554497-49673110', '3:49664859-49671529', '3:49671296-49671439'], ['dORF_56_XM_054328991.1', '-', '14:93237570-93333016', '14:93237569-93242252', '14:93237814-93237816'], ['dORF_2_NM_001321690.1', '+', '2:200811609-200827318', '2:200822198-200823831', '2:200822250-200822264'], ['dORF_16_NM_000719.7', '+', '12:1970800-2697930', '12:2691219-2697930', '12:2692748-2692837'], ['dORF_4_XM_005252591.4', '+', '10:18140444-18543537', '10:18539744-18543537', '10:18540167-18540178'], ['dORF_1_NM_172128.3', '-', '4:113451052-113761718', '4:113451051-113454491', '4:113454422-113454445'], ['dORF_5_XM_047426959.1', '+', '11:34051751-34102590', '11:34099387-34102590', '11:34099615-34099617'], ['dORF_2_XM_047426959.1', '+', '11:34051751-34102590', '11:34099387-34102590', '11:34099444-34099482'], ['dORF_70_XM_054379598.1', '-', '16:80597927-80805017', '16:80597926-80604367', '16:80597940-80597990'], ['dORF_7_NM_001300829.2', '+', '19:1269352-1274860', '19:1272463-1274790', '19:1273107-1273151'], ['dORF_10_XM_047447764.1', '-', '3:33496265-33718234', '3:33496264-33498610', '3:33497648-33497650'], ['dORF_7_NM_001204062.2', '+', '10:22316408-22331464', '10:22329562-22331464', '10:22330229-22330234'], ['dORF_11_XM_017018086.2', '-', '11:85658010-85682843', '11:85658009-85663565', '11:85661713-85661778'], ['dORF_12_XM_017018086.2', '-', '11:85658010-85682843', '11:85658009-85663565', '11:85661713-85661727'], ['dORF_16_XM_054359700.1', '-', '8:112222948-113436919', '8:112222947-112224750', '8:112223044-112223058'], ['dORF_6_NM_001288715.1', '-', '5:10971856-11904426', '5:10971859-10973432', '5:10972861-10972950'], ['dORF_17_NM_080759.6', '-', '13:71437986-71867184', '13:71437985-71440634', '13:71439331-71439339'], ['dORF_16_XM_054313463.1', '-', '16:57471942-57487302', '16:57471941-57473356', '16:57471985-57472002'], ['dORF_32_XM_024453007.2', '+', '2:26848015-26950331', '2:26947015-26950331', '2:26949906-26949911'], ['dORF_18_NM_001364157.2', '-', '5:158695940-159099896', '5:158695939-158699122', '5:158697027-158697071'], ['dORF_7_NM_001964.3', '+', '5:138465499-138469283', '5:138468101-138469283', '5:138468510-138468518'], ['dORF_45_NM_001243513.1', '-', '1:216503266-217137682', '1:216503265-216506918', '1:216503363-216503374'], ['dORF_22_NM_016605.3', '+', '5:138337577-138349709', '5:138346979-138349709', '5:138349588-138349668'], ['dORF_23_NM_203301.4', '-', '14:39397704-39432413', '14:39397703-39399495', '14:39397786-39397788'], ['dORF_14_NM_003868.3', '+', 'X:77447409-77457258', 'X:77456542-77457258', 'X:77457191-77457202'], ['dORF_30_NM_001346114.2', '-', '5:131641734-131796997', '5:131641733-131644664', '5:131641769-131641804'], ['dORF_3_XM_054339463.1', '-', '1:77944075-77979481', '1:77944074-77948745', '1:77948560-77948601'], ['dORF_6_NM_138426.4', '+', '7:7968816-8089060', '7:8086558-8089060', '7:8086899-8087018'], ['dORF_9_NM_022130.4', '-', '5:32124736-32174299', '5:32124735-32126191', '5:32124752-32124898'], ['dORF_2_NM_007325.5', '+', 'X:123184298-123490895', 'X:123488731-123490895', 'X:123488823-123488825'], ['dORF_2_NM_024503.5', '-', '1:41506385-42035914', '1:41506384-41510430', '1:41510281-41510283'], ['dORF_11_XM_054358064.1', '-', '7:26189947-26200726', '7:26189946-26192495', '7:26191092-26191100'], ['dORF_17_NM_001363050.1', '-', '8:74234720-74321520', '8:74234719-74236591', '8:74234863-74234931'], ['dORF_20_NM_001363050.1', '-', '8:74234720-74321520', '8:74234719-74236591', '8:74234857-74234862'], ['dORF_19_NM_002253.4', '-', '4:55078501-55125575', '4:55078500-55079920', '4:55078567-55078569'], ['dORF_13_XM_011532501.3', '+', '2:23385199-23708586', '2:23706684-23708586', '2:23708250-23708342'], ['dORF_2_XM_047430302.1', '+', '13:113297259-113323652', '13:113322441-113323652', '13:113322609-113322620'], ['dORF_54_XM_054331699.1', '-', '12:12116045-12267024', '12:12116044-12121105', '12:12116226-12116309'], ['dORF_11_NM_004631.5', '-', '1:53242384-53328050', '1:53242383-53246997', '1:53246779-53246790'], ['dORF_4_NM_004631.5', '-', '1:53242384-53328050', '1:53242383-53246997', '1:53246785-53246925'], ['dORF_3_NM_001350216.3', '+', '1:69567942-70144344', '1:70121796-70144344', '1:70121948-70121989'], ['dORF_5_NM_015578.4', '+', '19:34172524-34229268', '19:34227408-34229268', '19:34227931-34227969'], ['dORF_8_NM_001351625.3', '-', '4:86010425-86594054', '4:86010424-86017338', '4:86016810-86016857'], ['dORF_8_XM_054346559.1', '+', '3:152243652-152465760', '3:152462404-152465760', '3:152462961-152462969'], ['dORF_8_NM_007358.4', '+', '1:93079303-93139056', '1:93137047-93139056', '1:93137685-93137693'], ['dORF_9_NM_007358.4', '+', '1:93079303-93139056', '1:93137047-93139056', '1:93137691-93137693'], ['dORF_45_NM_001329851.3', '-', '2:1789133-2331255', '2:1789155-1792288', '2:1789171-1789230'], ['dORF_7_NM_021076.4', '+', '22:29480238-29491370', '22:29490723-29491370', '22:29491075-29491080'], ['dORF_7_NM_001135659.3', '-', '2:49918523-51032112', '2:49918522-49921923', '2:49921178-49921228'], ['dORF_2_XM_047418839.1', '-', '6:32184753-32190182', '6:32184752-32186361', '6:32185749-32185889'], ['dORF_12_XM_047430276.1', '+', '13:57630128-57729291', '13:57725314-57729291', '13:57726456-57726689'], ['dORF_14_XM_047430276.1', '+', '13:57630128-57729291', '13:57725314-57729291', '13:57726523-57726528'], ['dORF_4_XM_011535918.4', '-', '6:78934439-79078234', '6:78934438-78940672', '6:78940036-78940113'], ['dORF_22_NM_001382323.2', '+', '11:125164771-125433369', '11:125431412-125433369', '11:125433165-125433326'], ['dORF_23_NM_001382323.2', '+', '11:125164771-125433369', '11:125431412-125433369', '11:125433321-125433326'], ['dORF_1_XM_054331366.1', '-', '6:30600433-30617223', '6:30600432-30601528', '6:30601392-30601463'], ['dORF_2_XM_017001874.2', '+', '1:116909936-116990333', '1:116986987-116990333', '1:116987154-116987162'], ['dORF_32_NM_130393.3', '-', '9:8314266-10612982', '9:8314265-8317853', '9:8314827-8314832'], ['dORF_19_NM_001020658.2', '-', '1:30931526-31065697', '1:30931525-30933190', '1:30931576-30931584'], ['dORF_31_NM_206854.3', '+', '6:163414738-163578572', '6:163564730-163578572', '6:163566694-163566702'], ['dORF_5_NM_001378107.1', '+', '2:135531504-135725249', '2:135724312-135725249', '2:135724558-135724560'], ['dORF_26_NM_001100588.3', '-', '9:122844576-122905339', '9:122844575-122849606', '9:122847438-122847446'], ['dORF_22_XM_054361422.1', '-', '8:91954987-92103365', '8:91954986-91960221', '8:91958746-91958824'], ['dORF_23_XM_054361422.1', '-', '8:91954987-92103365', '8:91954986-91960221', '8:91958746-91958812'], ['dORF_19_NM_175634.3', '-', '8:91954987-92103365', '8:91954991-91960221', '8:91958826-91959068'], ['dORF_23_NM_175634.3', '-', '8:91954987-92103365', '8:91954991-91960221', '8:91958702-91958812'], ['dORF_3_NM_175634.3', '-', '8:91954987-92103365', '8:91954991-91960221', '8:91959829-91959954'], ['dORF_5_NM_175634.3', '-', '8:91954987-92103365', '8:91954991-91960221', '8:91959878-91959883'], ['dORF_1_XM_054375659.1', '+', '14:54565336-54793295', '14:54788964-54793295', '14:54788973-54789068'], ['dORF_9_NM_002971.6', '-', '3:18345397-18445572', '3:18345396-18349149', '3:18348634-18348645'], ['dORF_20_NM_001172509.2', '-', '2:199269520-199471246', '2:199269524-199272190', '2:199269786-199269788'], ['dORF_29_NM_001358351.3', '+', '15:47184109-47774208', '15:47771805-47774208', '15:47774183-47774185'], ['dORF_11_NM_015559.3', '+', '18:44680093-45068490', '18:45063718-45068490', '18:45065071-45065178'], ['dORF_7_NM_001372044.2', '+', '22:50672843-50733192', '22:50731332-50733192', '22:50732685-50732759'], ['dORF_5_XM_047426145.1', '-', '1:232397985-232630476', '1:232397984-232399106', '1:232398500-232398697'], ['dORF_6_XM_047426145.1', '-', '1:232397985-232630476', '1:232397984-232399106', '1:232398500-232398667'], ['dORF_1_NM_005413.4', '+', '2:44941722-44946051', '2:44944780-44946051', '2:44945169-44945180'], ['dORF_21_NM_001304421.2', '+', '3:57756329-57929993', '3:57927408-57929993', '3:57928966-57928989'], ['dORF_2_XM_054321826.1', '+', '19:10961050-11062253', '19:11061836-11062253', '19:11062236-11062238'], ['dORF_1_XM_054365273.1', '+', '10:104641310-105265222', '10:105263394-105265222', '10:105263441-105263521'], ['dORF_5_NM_001172712.1', '-', '2:173900795-173965682', '2:173906478-173909920', '2:173909734-173909736'], ['dORF_6_NM_001309444.2', '-', '5:151661116-151686895', '5:151661115-151663434', '5:151662548-151662628'], ['dORF_4_XM_047436525.1', '-', '17:57989058-58007226', '17:57989057-58005385', '17:58005038-58005097'], ['dORF_11_NM_001169117.2', '+', '4:26860861-27025361', '4:27022538-27025361', '4:27023324-27023326'], ['dORF_4_NM_182910.2', '+', '14:63761616-64226429', '14:64225546-64226473', '14:64226067-64226147'], ['dORF_1_XM_054356666.1', '+', '6:33418187-33453669', '6:33451926-33453669', '6:33452396-33452602'], ['dORF_63_XM_017024430.3', '+', '17:62966255-63427683', '17:63421975-63427683', '17:63427623-63427664'], ['dORF_1_XM_054325883.1', '-', '22:42160033-42343517', '22:42160032-42168632', '22:42161124-42168610'], ['dORF_2_XM_054325883.1', '-', '22:42160033-42343517', '22:42160032-42168632', '22:42161124-42161324'], ['dORF_1_XM_054325887.1', '-', '22:42160033-42343517', '22:42160032-42161320', '22:42161124-42161300'], ['dORF_26_NM_001098816.3', '-', '11:78652849-79441010', '11:78652848-78658037', '11:78656314-78656349'], ['dORF_6_XM_047436271.1', '+', '1:36207657-36305337', '1:36304037-36305337', '1:36304634-36304714'], ['dORF_13_XM_011518972.4', '+', '9:79571985-79726862', '9:79725164-79726862', '9:79726261-79726290'], ['dORF_20_NM_001195071.1', '-', '22:38216415-38272990', '22:38219310-38221448', '22:38219377-38219394'], ['dORF_17_XM_054316816.1', '+', '17:77957577-78108815', '17:78104865-78108815', '17:78106727-78106735'], ['dORF_14_NM_001080430.4', '-', '16:52436436-52547782', '16:52436436-52439204', '16:52438093-52438104'], ['dORF_23_NM_014112.5', '-', '8:115408516-115668955', '8:115408515-115414002', '8:115412444-115412446'], ['dORF_21_XM_054347743.1', '-', '3:33388356-33441379', '3:33388355-33390310', '3:33388539-33388661'], ['dORF_9_XM_047416346.1', '-', '4:95162524-95548953', '4:95162523-95169213', '4:95168817-95168852'], ['dORF_3_XM_047419923.1', '+', '7:65865792-65959538', '7:65954433-65959538', '7:65954650-65954658'], ['dORF_19_NM_017590.6', '+', '22:41301545-41360127', '22:41357449-41360127', '22:41359341-41359346'], ['dORF_17_NM_001278244.1', '-', '9:72351433-72365188', '9:72351444-72355932', '9:72354859-72354864'], ['dORF_38_XM_047434169.1', '-', '16:72782905-73891910', '16:72782904-72787143', '16:72783114-72783191'], ['dORF_28_XM_017027246.3', '+', '19:36214743-36238748', '19:36237339-36238748', '19:36238715-36238735'], ['dORF_23_XM_011526270.4', '-', '18:76357702-76496399', '18:76357701-76362477', '18:76361102-76361131'], ['dORF_75_XM_011521543.4', '+', '15:89998327-90082171', '15:90074670-90082171', '15:90082093-90082161']])

#create loop to take each SRA data in list in turn and run analysis in the loop
for i in range(len(SRA_Datasets)):
    #Create new temporary dataframe to hold results of each loop
    tempResultsDF = pd.DataFrame(columns=columns)
    #set dataset variable as dataset in loop
    dataset = SRA_Datasets[i]
    #Add dataset to the temp results df
    tempResultsDF['Dataset'] = [dataset]
    #uncompress sorted bam file for subsequent analysis
    os.system("gzip -d ./hisat2_alignments/" + dataset + "_sorted.bam")
    #Filter BAM file alignments based on the regions they are aligned to
    #loop through each of the regions in list for each aligned dataset
    for x in range(len(Housekeeping_Region_list)):
        #take details from housekeeping region list ready to use
        gene = Housekeeping_Region_list[x][0]
        strand = Housekeeping_Region_list[x][1]
        gene_region = Housekeeping_Region_list[x][2]
        three_prime_region = Housekeeping_Region_list[x][3]
        #filter BAM file to only include alignments within specified region
        samtools_view_gene = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + gene_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_gene_HAS.bam")
        os.system(samtools_view_gene)
        samtools_view_three = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + three_prime_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_three_prime_HAS.bam")
        os.system(samtools_view_three)
        
        #Generate results file with number of reads mapped to regions of interest
        bam_stats_gene = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_gene_HAS.bam > ./BAM_stats_results/" + dataset + "_" + gene + "_gene_HAS.txt")
        os.system(bam_stats_gene)
        bam_stats_three = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_three_prime_HAS.bam > ./BAM_stats_results/" + dataset + "_" + gene + "_three_prime_HAS.txt")
        os.system(bam_stats_three)
        
        #open BAM results for RP alignments to gene in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + gene + "_gene_HAS.txt")
        with open(open_string) as BAM_results:
            Temp_BAM_results = BAM_results.read()
        #Go through results line by line    
        for lines in Temp_BAM_results.split("\n"):
            #if the region of interest is on the +ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            if '+' in strand:
                if "Reads map to '+':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
            #if the region of interest is on the -ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            else:
                if "Reads map to '-':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
        
        #use regular expression to take number from line string variable and convert to integer
        Temp_reads = re.search("[0-9]{1,}", Temp_reads)
        Temp_reads = Temp_reads.group(0)
        Temp_reads = int(Temp_reads)
        
        #make sure results go into temp df column corresponding to gene in loop
        column_title = str(gene + " RP Alignments to Gene")
        tempResultsDF[column_title] = [Temp_reads]
        
        #remove temporary variable
        del(Temp_reads)
        
        #open BAM results for RP alignments to 3UTR in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + gene + "_three_prime_HAS.txt")
        with open(open_string) as BAM_results:
            Temp_BAM_results = BAM_results.read()
        #Go through results line by line    
        for lines in Temp_BAM_results.split("\n"):
            #if the region of interest is on the +ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            if '+' in strand:
                if "Reads map to '+':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
            #if the region of interest is on the -ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            else:
                if "Reads map to '-':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
        
        #use regular expression to take number from line string variable and convert to integer
        Temp_reads = re.search("[0-9]{1,}", Temp_reads)
        Temp_reads = Temp_reads.group(0)
        Temp_reads = int(Temp_reads)
        
        #make sure results go into temp df column corresponding to gene in loop
        column_title = str(gene + " RP Alignments to 3UTR")
        tempResultsDF[column_title] = [Temp_reads]
        
        #remove temporary variable
        del(Temp_reads)
        
    #loop through each of the regions in list for each aligned dataset
    for x in range(len(dORF_Region_List_HAS)):
        #take details from housekeeping region list ready to use
        dORF = dORF_Region_List_HAS[x][0]
        strand = dORF_Region_List_HAS[x][1]
        gene_region = dORF_Region_List_HAS[x][2]
        three_prime_region = dORF_Region_List_HAS[x][3]
        dORF_region = dORF_Region_List_HAS[x][4]
        #filter BAM file to only include alignments withn specified region
        samtools_view_gene = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + gene_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_gene_HAS.bam")
        os.system(samtools_view_gene)
        samtools_view_three = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + three_prime_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_three_prime_HAS.bam")
        os.system(samtools_view_three)
        samtools_view_dORF = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + dORF_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_dORF_HAS.bam")
        os.system(samtools_view_dORF)
        
        #Generate results file with number of reads mapped to regions of interest
        bam_stats_gene = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_gene_HAS.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_gene_HAS.txt")
        os.system(bam_stats_gene)
        bam_stats_three = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_three_prime_HAS.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_three_prime_HAS.txt")
        os.system(bam_stats_three)
        bam_stats_dORF = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_dORF_HAS.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_dORF_HAS.txt")
        os.system(bam_stats_dORF)
        
        #open BAM results for RP alignments to gene in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_gene_HAS.txt")
        with open(open_string) as BAM_results:
            Temp_BAM_results = BAM_results.read()
        #Go through results line by line    
        for lines in Temp_BAM_results.split("\n"):
            #if the region of interest is on the +ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            if '+' in strand:
                if "Reads map to '+':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
            #if the region of interest is on the -ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            else:
                if "Reads map to '-':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
        
        #use regular expression to take number from line string variable and convert to integer
        Temp_reads = re.search("[0-9]{1,}", Temp_reads)
        Temp_reads = Temp_reads.group(0)
        Temp_reads = int(Temp_reads)
        
        #make sure results go into temp df column corresponding to gene in loop
        column_title = str(dORF + " RP Alignments to Gene")
        tempResultsDF[column_title] = [Temp_reads]
        
        #remove temporary variable
        del(Temp_reads)
        
        #open BAM results for RP alignments to 3UTR in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_three_prime_HAS.txt")
        with open(open_string) as BAM_results:
            Temp_BAM_results = BAM_results.read()
        #Go through results line by line    
        for lines in Temp_BAM_results.split("\n"):
            #if the region of interest is on the +ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            if '+' in strand:
                if "Reads map to '+':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
            #if the region of interest is on the -ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            else:
                if "Reads map to '-':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
        
        #use regular expression to take number from line string variable and convert to integer
        Temp_reads = re.search("[0-9]{1,}", Temp_reads)
        Temp_reads = Temp_reads.group(0)
        Temp_reads = int(Temp_reads)
        
        #make sure results go into temp df column corresponding to gene in loop
        column_title = str(dORF + " RP Alignments to 3UTR")
        tempResultsDF[column_title] = [Temp_reads]
        
        #remove temporary variable
        del(Temp_reads)

        #open BAM results for RP alignments to 3UTR in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_dORF_HAS.txt")
        with open(open_string) as BAM_results:
            Temp_BAM_results = BAM_results.read()
        #Go through results line by line    
        for lines in Temp_BAM_results.split("\n"):
            #if the region of interest is on the +ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            if '+' in strand:
                if "Reads map to '+':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
            #if the region of interest is on the -ive strand then pull out the line matching reads aligned to this strand and take line as string variable
            else:
                if "Reads map to '-':" in lines:
                    Temp_reads = lines.strip()
                else:
                    continue
        
        #use regular expression to take number from line string variable and convert to integer
        Temp_reads = re.search("[0-9]{1,}", Temp_reads)
        Temp_reads = Temp_reads.group(0)
        Temp_reads = int(Temp_reads)
        
        #make sure results go into temp df column corresponding to gene in loop
        column_title = str(dORF + " RP Alignments to dORF")
        tempResultsDF[column_title] = [Temp_reads]
        
        #remove temporary variable
        del(Temp_reads)
        
    #show filtering of bam and reporting reads aligned has been completed
    print(dataset, ": BAM file filtered and results reported")
    #compress sorted bam file to save storage
    os.system("gzip ./hisat2_alignments/" + dataset + "_sorted.bam")
    
    #write results of temporary results dataframe into created .csv results file
    tempResultsDF.to_csv('Part6_RP_alignment_gene_3UTR_HAS_dORF_results.csv', header=None, mode='a')
    print(dataset, ": Results added to final dataframe")
        
        
