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
columns = ['Dataset', 'ACTB RP Alignments to Gene', 'ACTB RP Alignments to 3UTR', 'TUBB RP Alignments to Gene', 'TUBB RP Alignments to 3UTR', 'GAPDH RP Alignments to Gene', 'GAPDH RP Alignments to 3UTR', '3_HOXA11-201_ENSG00000005073_ENST00000006015 RP Alignments to Gene', '3_HOXA11-201_ENSG00000005073_ENST00000006015 RP Alignments to 3UTR', '3_HOXA11-201_ENSG00000005073_ENST00000006015 RP Alignments to dORF', '22_PTPN18-201_ENSG00000072135_ENST00000175756 RP Alignments to Gene', '22_PTPN18-201_ENSG00000072135_ENST00000175756 RP Alignments to 3UTR', '22_PTPN18-201_ENSG00000072135_ENST00000175756 RP Alignments to dORF', '142_ACYP1-201_ENSG00000119640_ENST00000238618 RP Alignments to Gene', '142_ACYP1-201_ENSG00000119640_ENST00000238618 RP Alignments to 3UTR', '142_ACYP1-201_ENSG00000119640_ENST00000238618 RP Alignments to dORF', '226_BTG1-201_ENSG00000133639_ENST00000256015 RP Alignments to Gene', '226_BTG1-201_ENSG00000133639_ENST00000256015 RP Alignments to 3UTR', '226_BTG1-201_ENSG00000133639_ENST00000256015 RP Alignments to dORF', '548_KRT86-201_ENSG00000170442_ENST00000293525 RP Alignments to Gene', '548_KRT86-201_ENSG00000170442_ENST00000293525 RP Alignments to 3UTR', '548_KRT86-201_ENSG00000170442_ENST00000293525 RP Alignments to dORF', '584_LSM6-201_ENSG00000164167_ENST00000296581 RP Alignments to Gene', '584_LSM6-201_ENSG00000164167_ENST00000296581 RP Alignments to 3UTR', '584_LSM6-201_ENSG00000164167_ENST00000296581 RP Alignments to dORF', '846_SPINK6-201_ENSG00000178172_ENST00000325630 RP Alignments to Gene', '846_SPINK6-201_ENSG00000178172_ENST00000325630 RP Alignments to 3UTR', '846_SPINK6-201_ENSG00000178172_ENST00000325630 RP Alignments to dORF', '849_OAZ2-201_ENSG00000180304_ENST00000326005 RP Alignments to Gene', '849_OAZ2-201_ENSG00000180304_ENST00000326005 RP Alignments to 3UTR', '849_OAZ2-201_ENSG00000180304_ENST00000326005 RP Alignments to dORF', '921_TNFAIP2-201_ENSG00000185215_ENST00000333007 RP Alignments to Gene', '921_TNFAIP2-201_ENSG00000185215_ENST00000333007 RP Alignments to 3UTR', '921_TNFAIP2-201_ENSG00000185215_ENST00000333007 RP Alignments to dORF', '1275_MOSPD1-203_ENSG00000101928_ENST00000370783 RP Alignments to Gene', '1275_MOSPD1-203_ENSG00000101928_ENST00000370783 RP Alignments to 3UTR', '1275_MOSPD1-203_ENSG00000101928_ENST00000370783 RP Alignments to dORF', '1400_NINJ1-201_ENSG00000131669_ENST00000375446 RP Alignments to Gene', '1400_NINJ1-201_ENSG00000131669_ENST00000375446 RP Alignments to 3UTR', '1400_NINJ1-201_ENSG00000131669_ENST00000375446 RP Alignments to dORF', '1408_IRS2-201_ENSG00000185950_ENST00000375856 RP Alignments to Gene', '1408_IRS2-201_ENSG00000185950_ENST00000375856 RP Alignments to 3UTR', '1408_IRS2-201_ENSG00000185950_ENST00000375856 RP Alignments to dORF', '1548_HAP1-204_ENSG00000173805_ENST00000393939 RP Alignments to Gene', '1548_HAP1-204_ENSG00000173805_ENST00000393939 RP Alignments to 3UTR', '1548_HAP1-204_ENSG00000173805_ENST00000393939 RP Alignments to dORF', '1575_SLC43A3-202_ENSG00000134802_ENST00000395123 RP Alignments to Gene', '1575_SLC43A3-202_ENSG00000134802_ENST00000395123 RP Alignments to 3UTR', '1575_SLC43A3-202_ENSG00000134802_ENST00000395123 RP Alignments to dORF', '1663_SZRD1-202_ENSG00000055070_ENST00000401089 RP Alignments to Gene', '1663_SZRD1-202_ENSG00000055070_ENST00000401089 RP Alignments to 3UTR', '1663_SZRD1-202_ENSG00000055070_ENST00000401089 RP Alignments to dORF', '1711_SNRPE-202_ENSG00000182004_ENST00000414487 RP Alignments to Gene', '1711_SNRPE-202_ENSG00000182004_ENST00000414487 RP Alignments to 3UTR', '1711_SNRPE-202_ENSG00000182004_ENST00000414487 RP Alignments to dORF', '1770_TAPBP-234_ENSG00000231925_ENST00000434618 RP Alignments to Gene', '1770_TAPBP-234_ENSG00000231925_ENST00000434618 RP Alignments to 3UTR', '1770_TAPBP-234_ENSG00000231925_ENST00000434618 RP Alignments to dORF', '1833_C1orf52-203_ENSG00000162642_ENST00000471115 RP Alignments to Gene', '1833_C1orf52-203_ENSG00000162642_ENST00000471115 RP Alignments to 3UTR', '1833_C1orf52-203_ENSG00000162642_ENST00000471115 RP Alignments to dORF', '1849_POLE4-205_ENSG00000115350_ENST00000483063 RP Alignments to Gene', '1849_POLE4-205_ENSG00000115350_ENST00000483063 RP Alignments to 3UTR', '1849_POLE4-205_ENSG00000115350_ENST00000483063 RP Alignments to dORF', '1975_LIN52-205_ENSG00000205659_ENST00000555028 RP Alignments to Gene', '1975_LIN52-205_ENSG00000205659_ENST00000555028 RP Alignments to 3UTR', '1975_LIN52-205_ENSG00000205659_ENST00000555028 RP Alignments to dORF', '2021_COX6B2-205_ENSG00000160471_ENST00000588572 RP Alignments to Gene', '2021_COX6B2-205_ENSG00000160471_ENST00000588572 RP Alignments to 3UTR', '2021_COX6B2-205_ENSG00000160471_ENST00000588572 RP Alignments to dORF', '2024_NFIC-206_ENSG00000141905_ENST00000589123 RP Alignments to Gene', '2024_NFIC-206_ENSG00000141905_ENST00000589123 RP Alignments to 3UTR', '2024_NFIC-206_ENSG00000141905_ENST00000589123 RP Alignments to dORF', '2043_C19orf53-206_ENSG00000104979_ENST00000593274 RP Alignments to Gene', '2043_C19orf53-206_ENSG00000104979_ENST00000593274 RP Alignments to 3UTR', '2043_C19orf53-206_ENSG00000104979_ENST00000593274 RP Alignments to dORF', '2082_MARCKS-201_ENSG00000277443_ENST00000612661 RP Alignments to Gene', '2082_MARCKS-201_ENSG00000277443_ENST00000612661 RP Alignments to 3UTR', '2082_MARCKS-201_ENSG00000277443_ENST00000612661 RP Alignments to dORF', '103_SUPT4H1-201_ENSG00000213246_ENST00000225504 RP Alignments to Gene', '103_SUPT4H1-201_ENSG00000213246_ENST00000225504 RP Alignments to 3UTR', '103_SUPT4H1-201_ENSG00000213246_ENST00000225504 RP Alignments to dORF', '209_PXDN-201_ENSG00000130508_ENST00000252804 RP Alignments to Gene', '209_PXDN-201_ENSG00000130508_ENST00000252804 RP Alignments to 3UTR', '209_PXDN-201_ENSG00000130508_ENST00000252804 RP Alignments to dORF']
ResultsDF = pd.DataFrame(columns=columns)

#Save Results DF as empty .csv with headings for columns
ResultsDF.to_csv('Part6_RP_alignment_gene_3UTR_26MS_dORF_results.csv')
print('Empty Results Dataframe .csv created')

#list of regions to look at RP read alignment
Housekeeping_Region_list = list([['ACTB', '-', '7:5527168-5530581', '7:5527167-5527727'], ['TUBB', '+', '6:30720372-30725402', '6:30724417-30725402'], ['GAPDH', '+', '12:6534537-6538351', '12:6538190-6538351']])
dORF_Region_List_26MS = list([['3_HOXA11-201_ENSG00000005073_ENST00000006015', '-', '7:27181177-27185212', '7:27181534-27182775', '7:27182649-27182726'], ['22_PTPN18-201_ENSG00000072135_ENST00000175756', '+', '2:130356076-130375385', '2:130373244-130375389', '2:130373735-130373860'], ['142_ACYP1-201_ENSG00000119640_ENST00000238618', '-', '14:75053263-75064004', '14:75053240-75053423', '14:75053244-75053300'], ['226_BTG1-201_ENSG00000133639_ENST00000256015', '-', '12:92140298-92145826', '12:92140297-92144059', '12:92143962-92144045'], ['548_KRT86-201_ENSG00000170442_ENST00000293525', '+', '12:52274665-52309143', '12:52308605-52309143', '12:52309000-52309065'], ['584_LSM6-201_ENSG00000164167_ENST00000296581', '+', '4:146175738-146191515', '4:146189676-146191515', '4:146189814-146189873'], ['846_SPINK6-201_ENSG00000178172_ENST00000325630', '+', '5:148202798-148215117', '5:148214970-148215117', '5:148214980-148215054'], ['849_OAZ2-201_ENSG00000180304_ENST00000326005', '-', '15:64687594-64703261', '15:64687592-64688683', '15:64688510-64688671'], ['921_TNFAIP2-201_ENSG00000185215_ENST00000333007', '+', '14:103121489-103137419', '14:103135380-103137419', '14:103135741-103135983'], ['1275_MOSPD1-203_ENSG00000101928_ENST00000370783', '-', 'X:134887652-134915237', 'X:134887645-134889140', 'X:134888605-134888688'], ['1400_NINJ1-201_ENSG00000131669_ENST00000375446', '-', '9:93121516-93134231', '9:93121508-93122507', '9:93121895-93122128'], ['1408_IRS2-201_ENSG00000185950_ENST00000375856', '-', '13:109752715-109786563', '13:109752717-109756283', '13:109755764-109755997'], ['1548_HAP1-204_ENSG00000173805_ENST00000393939', '-', '17:41717759-41734626', '17:41720832-41724680', '17:41721488-41721547'], ['1575_SLC43A3-202_ENSG00000134802_ENST00000395123', '-', '11:57406974-57427560', '11:57406973-57407771', '11:57407432-57407557'], ['1663_SZRD1-202_ENSG00000055070_ENST00000401089', '+', '1:16367262-16398125', '1:16395160-16398125', '1:16395160-16395294'], ['1711_SNRPE-202_ENSG00000182004_ENST00000414487', '+', '1:203861619-203871132', '1:203869952-203870530', '1:203870066-203870143'], ['1770_TAPBP-234_ENSG00000231925_ENST00000434618', '-', '6:33299714-33314058', '6:33299713-33301739', '6:33299742-33299792'], ['1833_C1orf52-203_ENSG00000162642_ENST00000471115', '-', '1:85249973-85259642', '1:85251987-85252608', '1:85252447-85252602'], ['1849_POLE4-205_ENSG00000115350_ENST00000483063', '+', '2:74958663-74970108', '2:74969442-74970108', '2:74969653-74969703'], ['1975_LIN52-205_ENSG00000205659_ENST00000555028', '+', '14:74084976-74201473', '14:74198997-74201215', '14:74199015-74199197'], ['2021_COX6B2-205_ENSG00000160471_ENST00000588572', '-', '19:55349724-55354699', '19:55349325-55353700', '19:55350577-55350783'], ['2024_NFIC-206_ENSG00000141905_ENST00000589123', '+', '19:3359650-3469197', '19:3462789-3469197', '19:3463624-3463692'], ['2043_C19orf53-206_ENSG00000104979_ENST00000593274', '+', '19:13774476-13778753', '19:13778218-13778442', '19:13778264-13778425'], ['2082_MARCKS-201_ENSG00000277443_ENST00000612661', '+', '6:113857365-113863455', '6:113860599-113863451', '6:113861644-113861712'], ['103_SUPT4H1-201_ENSG00000213246_ENST00000225504', '-', '17:58345198-58352181', '17:58345197-58346225', '17:58345745-58345996'], ['209_PXDN-201_ENSG00000130508_ENST00000252804', '-', '2:1631907-1744881', '2:1631906-1634183', '2:1634019-1634126']])

#create loop to take each SRA data in list in turn and run analysis in the loop
for i in range(len(SRA_Datasets)):
    #Create new temporary dataframe to hold results of each loop
    tempResultsDF = pd.DataFrame(columns=columns)
    #set dataset variable as dataset in loop
    dataset = SRA_Datasets[i]
    #Add dataset to the temp results df
    tempResultsDF['Dataset'] = [dataset]
    #extract SRA dataset
    extract_sra_dataset = str("fasterq-dump --outdir ./Extracted_SRA_Datasets " + dataset + " >/dev/null 2>&1")
    os.system(extract_sra_dataset)
    #show extraction complete
    print(dataset, ': Extraction Completed')
    
    #Run quality check of dataset with fastqc
    run_fastqc = str("../FastQC/fastqc --outdir ./Extracted_SRA_Datasets ./Extracted_SRA_Datasets/" + dataset + ".fastq >/dev/null 2>&1")
    os.system(run_fastqc)
    #unzip fastqc report
    unzip_fastqc = str("unzip -d ./Extracted_SRA_Datasets/ ./Extracted_SRA_Datasets/" + dataset + "_fastqc.zip >/dev/null 2>&1")
    os.system(unzip_fastqc)
    #remove files no longer needed
    remove_zip = str("rm ./Extracted_SRA_Datasets/" + dataset + "_fastqc.zip >/dev/null 2>&1")
    os.system(remove_zip)
    remove_html = str("rm ./Extracted_SRA_Datasets/" + dataset + "_fastqc.html >/dev/null 2>&1")
    os.system(remove_html)
    
    #Trim each different adapter sequences from fastq dataset with overlap at least 5 bases and Trim based on quality, n content and length of reads, also remove previous trim to save space
    trim_adapters_1 = str("../TrimGalore/trim_galore --illumina --stringency 5 --quality 30 --length 25 --trim-n --dont_gzip --basename " + dataset + "_i --output_dir ./Trimmed_SRA_Datasets ./Extracted_SRA_Datasets/" + dataset + ".fastq >/dev/null 2>&1")
    os.system(trim_adapters_1)
    trim_adapters_2 = str("../TrimGalore/trim_galore --nextera --stringency 5 --quality 30 --length 25 --trim-n --dont_gzip --basename " + dataset + "_in --output_dir ./Trimmed_SRA_Datasets ./Trimmed_SRA_Datasets/" + dataset + "_i_trimmed.fq >/dev/null 2>&1")
    os.system(trim_adapters_2)
    os.system("rm ./Trimmed_SRA_Datasets/" + dataset + "_i_trimmed.fq")
    trim_adapters_3 = str("../TrimGalore/trim_galore --small_rna --stringency 5 --quality 30 --length 25 --trim-n --max_n 0 --dont_gzip --basename " + dataset + "_ins --output_dir ./Trimmed_SRA_Datasets ./Trimmed_SRA_Datasets/" + dataset + "_in_trimmed_trimmed.fq >/dev/null 2>&1")
    os.system(trim_adapters_3)
    os.system("rm ./Trimmed_SRA_Datasets/" + dataset + "_in_trimmed_trimmed.fq")
    #Use hard trim if no adapters removed trims reads to 35 bases max length
    trim_length = str("../TrimGalore/trim_galore --dont_gzip --hardtrim5 35 --output_dir ./Trimmed_SRA_Datasets ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.fq >/dev/null 2>&1")
    os.system(trim_length)
    os.system("rm ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.fq")
    #compress the original extracted dataset to save space
    os.system("gzip ./Extracted_SRA_Datasets/" + dataset + ".fastq")
    #show trimming complete
    print(dataset, ': Trimming Completed')
    
    #Run quality check of trimmed dataset with fastqc
    run_fastqc_trim = str("../FastQC/fastqc --outdir ./Trimmed_SRA_Datasets ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime.fq >/dev/null 2>&1")
    os.system(run_fastqc_trim)
    #unzip fastqc report
    unzip_fastqc_trim = str("unzip -d ./Trimmed_SRA_Datasets ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime_fastqc.zip >/dev/null 2>&1")
    os.system(unzip_fastqc_trim)
    #remove files no longer needed
    remove_zip_trim = str("rm ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime_fastqc.zip >/dev/null 2>&1")
    os.system(remove_zip_trim)
    remove_html_trim = str("rm ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime_fastqc.html >/dev/null 2>&1")
    os.system(remove_html_trim)
    
    #Run hisat2 alignment between processed fastq dataset and human reference genome grch38 with high stringency
    hisat2 = str("../hisat2/hisat2 --n-ceil L,0.0,0.0 --score-min L,0.0,0.0 -x grch38 -U ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime.fq > ./hisat2_alignments/" + dataset + ".bam")
    os.system(hisat2)
    #compress trimmed fastq dataset to save space
    os.system("gzip ./Trimmed_SRA_Datasets/" + dataset + "_ins_trimmed_trimmed_trimmed.35bp_5prime.fq")
    #show hisat2 alignment complete
    print(dataset, ': GRCh38 Alignment Completed')
    
    #Sort and Index bam alignment file before it can be filtered then remove original BAM file
    sort = str("../samtools/samtools sort -u -o ./hisat2_alignments/" + dataset + "_sorted.bam ./hisat2_alignments/" + dataset + ".bam >/dev/null 2>&1")
    os.system(sort)
    index = str("../samtools/samtools index ./hisat2_alignments/" + dataset + "_sorted.bam >/dev/null 2>&1")
    os.system(index)    
    #save space by removing BAM file no longer needed
    remove_BAM = str("rm ./hisat2_alignments/" + dataset + ".bam >/dev/null 2>&1")
    os.system(remove_BAM)
    #Show bam file sorted and indexed
    print(dataset, ": BAM file sorted and indexed")
    
    #Filter BAM file alignments based on the regions they are aligned to
    #loop through each of the regions in list for each aligned dataset
    for x in range(len(Housekeeping_Region_list)):
        #take details from housekeeping region list ready to use
        gene = Housekeeping_Region_list[x][0]
        strand = Housekeeping_Region_list[x][1]
        gene_region = Housekeeping_Region_list[x][2]
        three_prime_region = Housekeeping_Region_list[x][3]
        #filter BAM file to only include alignments within specified region
        samtools_view_gene = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + gene_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_gene.bam")
        os.system(samtools_view_gene)
        samtools_view_three = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + three_prime_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_three_prime.bam")
        os.system(samtools_view_three)
        
        #Generate results file with number of reads mapped to regions of interest
        bam_stats_gene = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_gene.bam > ./BAM_stats_results/" + dataset + "_" + gene + "_gene.txt")
        os.system(bam_stats_gene)
        bam_stats_three = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + gene + "_three_prime.bam > ./BAM_stats_results/" + dataset + "_" + gene + "_three_prime.txt")
        os.system(bam_stats_three)
        
        #open BAM results for RP alignments to gene in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + gene + "_gene.txt")
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
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + gene + "_three_prime.txt")
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
    for x in range(len(dORF_Region_List_26MS)):
        #take details from housekeeping region list ready to use
        dORF = dORF_Region_List_26MS[x][0]
        strand = dORF_Region_List_26MS[x][1]
        gene_region = dORF_Region_List_26MS[x][2]
        three_prime_region = dORF_Region_List_26MS[x][3]
        dORF_region = dORF_Region_List_26MS[x][4]
        #filter BAM file to only include alignments withn specified region
        samtools_view_gene = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + gene_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_gene.bam")
        os.system(samtools_view_gene)
        samtools_view_three = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + three_prime_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_three_prime.bam")
        os.system(samtools_view_three)
        samtools_view_dORF = str("../samtools/samtools view -b -h ./hisat2_alignments/" + dataset + "_sorted.bam " + dORF_region + " > ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_dORF.bam")
        os.system(samtools_view_dORF)
        
        #Generate results file with number of reads mapped to regions of interest
        bam_stats_gene = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_gene.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_gene.txt")
        os.system(bam_stats_gene)
        bam_stats_three = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_three_prime.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_three_prime.txt")
        os.system(bam_stats_three)
        bam_stats_dORF = str("rseqc bam_stats -i ./filtered_hisat2_alignments/" + dataset + "_" + dORF + "_dORF.bam > ./BAM_stats_results/" + dataset + "_" + dORF + "_dORF.txt")
        os.system(bam_stats_dORF)
        
        #open BAM results for RP alignments to gene in readable format and convert to variable with STR of all text
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_gene.txt")
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
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_three_prime.txt")
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
        open_string = str("/users/payjt2/PhD/RP_dataset_dORF_alignment_analysis/BAM_stats_results/" + dataset + "_" + dORF + "_dORF.txt")
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
    tempResultsDF.to_csv('Part6_RP_alignment_gene_3UTR_26MS_dORF_results.csv', header=None, mode='a')
    print(dataset, ": Results added to final dataframe")
        
        
