#!/usr/bin/env python3
"""
workflow for mapping AmpliconSeq data and estimating DNA source including mosquito species and host
"""
#### Standard workflow imports and defs ###################################### 
import subprocess
import glob
import sys
import pandas as pd
import pandas.io.common
import gzip
import natsort as ns
import xlsxwriter
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', 200)
#### End standard workflow imports and defs ###################################


# SET THESE...
# import raw illumina read 1 files for each sample (This script assumes you have paired-end data and would need to be updated if you only have SE reads. 
R1s = glob.glob("aegy-LA-MD-*L001_R1_001.fastq.gz")  # add the full path if the files are not in the same directory
min_read_percent = 10 # percent threshold: blast worthy
min_identical_reads = 2 # percent threshold: blast worthy
read_len = 75
blast_output_format = "10"

########## Dictionaries ##############
# we do not have sufficient seq data to distinguish Canis species. I also create a shorter name for species.
name_conversion_dictionary = {"Culex quinquefasciatus mitochondrion":"quinq", "Culex pipiens pipiens mitochondrion":"pipiens","Aedes albopictus mitochondrion":"albo", "Aedes notoscriptus isolate MC04 mitochondrion":"noto", "Aedes aegypti strain LVP_AGWG mitochondrion":"aegypti", "Homo sapiens mitochondrion":"Human","Canis lupus lupus mitochondrion":"Canis","Canis lupus familiaris mitochondrion":"Canis", "Canis lupus mitochondrion":"Canis","Canis latrans mitochondrion":"Canis","Canis lupus chanco mitochondrion":"Canis","Canis lupus laniger mitochondrion":"Canis"}

# make genome_ID to common name dictionary
spp_D = {}
for line in open("/data/home/bmain/bin/BLAST_DB/mito/mito.nt"):
    if line[0]==">":
        i = line.strip().split(",")[0].split()
        ID = i[0].strip(">")
        spp = " ".join(i[1:])
        spp_D[ID] = spp

########## End of Dictionaries ##############


######## Start Functions #############
# Combine identical reads and filter > 10% (this is conservative) to remove error reads 
def group_reads(fq_read1, primer_seq,locus, size, threshold):
    df = pd.read_csv(fq_read1,engine='python', sep='\s',header = None,dtype = {0 :str})
    ids = list(df[df.index % 4 == 0][0])
    reads = list(df[df.index % 4 == 1][0])
    df1 = pd.DataFrame(list(zip(ids, reads)),columns =['id', 'read']) 
    fq_read2 = fq_read1.replace("R1","R2")
    df2 = pd.read_csv(fq_read2,engine='python', sep='\s',header = None,dtype = {0 :str})
    ids2 = list(df2[df2.index % 4 == 0][0])
    reads2 = list(df2[df2.index % 4 == 1][0])
    df2 = pd.DataFrame(list(zip(ids2, reads2)),columns =['id', 'read2']) 
    joined = pd.merge(df1, df2, on='id')
    locus_specific_DF = joined.loc[joined['read'].str.contains(primer_seq)]
    gap = size - (2*read_len) # read_len variable is define at top
    #locus_specific_DF['full_read'] = locus_specific_DF['read'] + gap*"N" + locus_specific_DF['read2']
    locus_specific_DF['full_read'] = locus_specific_DF['read'] + "NNNNNN" + locus_specific_DF['read2'] # UNLESS THE READS OVERLAP, Adding the proper size N gap does not appear to improved BLAST mapping. It breaks up the fragment anyway.
    locus_specific_DF = locus_specific_DF.groupby("full_read").size().sort_values(ascending=False).reset_index(name='count')
    Total = locus_specific_DF['count'].sum()
    locus_specific_DF['sum_reads'] = locus_specific_DF['count'].sum() 
    locus_specific_DF['percent'] = (100 * (locus_specific_DF['count'] / Total)).astype(int)
    locus_specific_DF = locus_specific_DF.loc[locus_specific_DF['count'] >= min_identical_reads] 
    locus_specific_DF['query_id'] = ">" + fq_read1.split("_")[0] + "_" + locus + "_" + locus_specific_DF['count'].astype(str) + "_" + Total.astype(str) + "_" + locus_specific_DF['percent'].astype(str) + "%" "__" + locus_specific_DF['full_read'].astype(str) + "_" 
    return(locus_specific_DF)

# blastn reads to my mitochondrial database. This have species known to be in CA and no prehistoric species. 
def blaster(read_DF, input_filename,amplicon_id):
    to_fasta_df = read_DF[['query_id','full_read']].stack()
    to_fasta_df = to_fasta_df.to_frame().reset_index()[[0]]
    to_fasta_df.to_csv("blaster_tmp.fa",index=False,header=False)
    blast = subprocess.call(["/data/home/bmain/bin/ncbi-blast-2.7.1+/bin/blastn", "-query", "/data/home/bmain/tarsalis/MiSeq_COAV_2019/marisa/blaster_tmp.fa", "-db", "/data/home/bmain/bin/BLAST_DB/mito/mito_removed_0322.nt", "-outfmt", blast_output_format, "-num_alignments", "8", "-out", "/data/home/bmain/tarsalis/MiSeq_COAV_2019/marisa/tmp_pandas.result", "-num_threads", "4"])
    try:
        results = pd.read_csv("tmp_pandas.result", header=None)
        results.columns = ["query_id", "matched_spp", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]
        results['bp_matching'] = results['alignment length'] - results['mismatches']
        results['hit'] = results['query_id'] + ":" + results['matched_spp']
        results = results.replace({"matched_spp": spp_D})
        results = results.groupby('hit').agg({"bp_matching":sum,"matched_spp":list}).reset_index()
        results['spp'] = results["matched_spp"].str[0] # in this case, the read is split, but both map to the same genome, soI arbitrarily take the first element.
        del results["matched_spp"]
        results = results.loc[results['bp_matching'] > 135]  # I determined this threshold by manually inspecting results from multiple samples. Below 135 there was less consensus.
        results = results[['hit','spp','bp_matching']]
        results.columns = ["query_id", "matched_spp", 'bp_matching']
    
    except pandas.io.common.EmptyDataError: 
        results = read_DF[['query_id']]
        results['matched_spp'] = "No_Hit"
        results['bp_matching'] = "No_Hit"
        results= results[["query_id", "matched_spp", 'bp_matching']]
        results['query_id'] = results['query_id'].str.replace(r'>', '')
    finally:
        results['amplicon_target'] = amplicon_id
        return(results)

def merge_DFs(short_DF, long_DF):
    long_DF["query_id"] = long_DF["query_id"].str.replace(r'>', '')
    no_genome_info = short_DF["query_id"].str.split(":", n=1,expand = True)
    short_DF["query_id"] = no_genome_info[0]
    plus_meta = pd.merge(short_DF, long_DF, on='query_id')
    return(plus_meta)


def ready_to_append(grouped_reads_DF,read1_file,target_locus):
    if len(grouped_reads_DF) > 0:
        blasted = blaster(grouped_reads_DF, read1_file, target_locus)
        if len(blasted) > 0:
            meta = merge_DFs(blasted, grouped_reads_DF)
            #return([blasted,meta])
            return(meta)
        else:
            sample_name = read1_file.split("_")[0] + "_" + target_locus 
            #d = {'query_id': [sample_name], 'matched_spp':["No_Hit"], 'bp_matching':["No_Hit"]}
            #no_hit_df = pd.DataFrame(data=d)
            d = {'query_id': [sample_name], 'matched_spp':["No_Hit"], 'bp_matching':["No_Hit"],'amplicon_target':["No_Hit"],'full_read':["No_Hit"], 'count':["No_Hit"], 'sum_reads':["No_Hit"], 'percent':["No_Hit"]}
            no_hit_df = pd.DataFrame(data=d)
            return(no_hit_df)
    else:
        sample_name = read1_file.split("_")[0] + "_" + target_locus 
        #d = {'query_id': [sample_name], 'matched_spp':["No_Hit"], 'bp_matching':["No_Hit"]}
        #no_hit_df = pd.DataFrame(data=d)
        d = {'query_id': [sample_name], 'matched_spp':["No_Hit"], 'bp_matching':["No_Hit"],'amplicon_target':["No_Hit"],'full_read':["No_Hit"], 'count':["No_Hit"], 'sum_reads':["No_Hit"], 'percent':["No_Hit"]}
        no_hit_df = pd.DataFrame(data=d)
        return(no_hit_df)

########## End function definition section ################
# Initialize dataframe
meta_df = pd.DataFrame(columns=['query_id', 'matched_spp', 'bp_matching', 'amplicon_target','full_read', 'count', 'sum_reads', 'percent'])

for R1 in R1s:
    R2 = R1.replace("R1","R2")
    # extract locus-specific reads for each F and R readset.
    mini_16s = "AAGACGAGAAGACCC"
    mini_16s_R = "TTGCGCTGTTATT"
    cytb = "CCATCCAACATCT"
    cytb_R = "CCTCGAATGAT"
    pgk2 = "AAGGTGAACGAAATGAT"
    pgk2 = "AAGGTGAACGAAATGAT"
    pgk2_R = "GTGTGATACCTG" 
    mtcox1 = "TGGGTCTCCTCC"
    mtcox1_R = "GAGCTCCTGATA"
    mtcox1_short = "CGATCAAGAGTAA"
    
    print("16s")
    df_16s = group_reads(R1,mini_16s, "16s",237, min_read_percent) # input expected amplicon size second to last
    add_16s = ready_to_append(df_16s, R1, "16s")
    meta_df = meta_df.append(add_16s)

    print("cytb")
    df_cytb = group_reads(R1,cytb, "cytb",357, min_read_percent) # input expected amplicon size second to last
    add_cytb = ready_to_append(df_cytb, R1, "cytb")
    meta_df = meta_df.append(add_cytb)
    

    print("mtcox1")
    df_mtcox1 = group_reads(R1,mtcox1, "mtcox1",423, min_read_percent) # input expected amplicon size second to last
    add_mtcox1 = ready_to_append(df_mtcox1, R1, "mtcox1")
    meta_df = meta_df.append(add_mtcox1)
    
    print("mtcox1_short")
    df_mtcox1_short = group_reads(R1,mtcox1_short, "mtcox1_short",178, min_read_percent) 
    add_mtcox1_short = ready_to_append(df_mtcox1_short, R1, "mtcox1_short")
    meta_df = meta_df.append(add_mtcox1_short)


meta_df.reset_index()
sample_name = meta_df['query_id'].str.split("_", n = 1, expand = True)
meta_df['sample'] = sample_name[0]
meta_df['sample'] = pd.Categorical(meta_df['sample'], ordered=True, categories= ns.natsorted(meta_df['sample'].unique()))
meta_df = meta_df.sort_values(['sample', 'amplicon_target', 'count'], ascending=[True,True,False])
meta_df = meta_df[['sample', 'amplicon_target', 'matched_spp', 'count', 'sum_reads', 'percent', 'bp_matching','full_read']]


# make easy to read calls sheet
meta_df.rename(columns={'sample':'sample_names'},inplace=True) # sample is a function
sample_list = meta_df.sample_names.unique()
target_list = list(meta_df.amplicon_target.unique())
target_list.remove("No_Hit")
meta_df = meta_df.replace({"matched_spp": name_conversion_dictionary})

calls={}
for s in sample_list:
    if s not in calls:
        calls[s] = {}
    for t in target_list:
        if t not in calls[s]:
            calls[s][t]=[]
        target = meta_df.loc[(meta_df['sample_names'] == s) & (meta_df['amplicon_target'] == t)]
        target['percent'] = target['percent'].astype(int)
        target = target.loc[target['percent'] >= min_read_percent]
        called_spp = list(target.matched_spp.unique())
        calls[s][t].extend(called_spp)

results = pd.DataFrame(data=calls)
results = results.T.reset_index()
for amplicon_target in target_list:
    results[amplicon_target]= results[amplicon_target].str.join("; ")

results.rename(columns={'index':'sample_ID'},inplace=True)


# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('marisa_version7.xlsx', engine='xlsxwriter')

# Write each dataframe to a different worksheet.
results.to_excel(writer, sheet_name='final_calls', index=False)
meta_df.to_excel(writer, sheet_name='meta_data', index=False)

# Close the Pandas Excel writer and output the Excel file.
writer.save()




