import subprocess
import glob
import sys

# import raw illumina read files for each sample 
R1s = glob.glob("*L001_R1_001.fastq.gz")
R1s.sort()

# This function identifies all unique sequences and counts abundance
def makeD(reads):
	haplotypes = {}
	counter=0
	for r in reads.split():
		r = r.strip('AAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG') # remove common bad sequence
		if r not in haplotypes:
			haplotypes[r]=0
		haplotypes[r]+=1
		counter+=1
	return [haplotypes, counter]

# This function formats each unique sequence to fasta and BLASTs the sequence to a modified (including only CA species + livestock and human/dog) mitochondrial BLAST database.
def blast(hap_dict, total_reads,filename,locus):	
	for hap in hap_dict:
		# I only BLAST sequences that are present at > 5% and have at least 3 copies
		if float(hap_dict[hap])/total_reads > 0.05 and hap_dict[hap]>2:
			
			# make fasta formated input file
			query = open("/data/home/bmain/tarsalis/multiplex/tmp.fa", "w")
			query.write(">%s" % File)
			query.write("\n")
			query.write(hap)
			query.close()
			# perform BLAST with 8 outputs and 40 threads.
			blast = subprocess.call(["/data/home/bmain/bin/ncbi-blast-2.7.1+/bin/blastn", "-query", "/data/home/bmain/tarsalis/multiplex/tmp.fa", "-db", "/data/home/bmain/bin/BLAST_DB/mito/custom_mito4_noRare.nt", "-outfmt", "3", "-num_alignments", "8", "-num_descriptions", "8", "-out", "/data/home/bmain/tarsalis/multiplex/tmp.result", "-num_threads", "40"])
			
			#open BLAST output
			f=open('/data/home/bmain/tarsalis/multiplex/tmp.result')
			lines=f.readlines()
			
			if len(lines)<5:
				print 'bad haplotype:', hap
				continue

			if "No" in lines[19].split(): # If no hits were found, the output says: "No results"
				continue #delete this if you want to see "no hit" info
				print filename, locus, hap,'hits:',hap_dict[hap], "_".join(lines[19].strip().split())
			else:
				Ps = [] # Pvalues
				spp_D={}
				for hit in range(20,25):
					if len(lines[hit].split()) < 2:
						break	# In cases where less than 8 BLAST results were identified, we stop extracting data
					
					exponent = int(lines[hit].split()[-1].split('e')[1]) # Make a list of the significance values that we use to sort the results with
					Ps.append(exponent)
				
				tophit= str(min(Ps)) # find the most significant Evalue
				spp_list = []
				spp_hits = []
				for hit in range(20,25):
					if len(lines[hit].split()) < 2:
						break
					x = lines[hit].split()	
					p = x[-1].split('e')[1] # Pvalue (Evalue)
					if int(p)<-60:		#imposed somewhat arbitrary Evalue threshold.
						if p == tophit:	# To narrow down the results, I only take the best hit/s
							spp_name = "_".join([x[1],x[2]])
							spp_list.append(spp_name)
							spp_list.append(p)
							spp_hits.append(spp_name)
				###print filename, locus, hap, 'hits:',hap_dict[hap], " ".join(spp_list)   #uncomment this if you want to see the detailed blast hit info
				return(spp_hits)			


for File in R1s:
	D={}
	D2=[]
	# extract locus-specific reads for each F and R readset.
	reads_16s = subprocess.check_output(["zgrep", "CGCCTGTTTATCAAAAACAT", File]) # 16S
	reads_cytB = subprocess.check_output(["zgrep", "CCATCCAACATCTCAGCATGATGAAA", File]) # cytB
	reads_16s_R2 = subprocess.check_output(["zgrep", "CGGTCTGAACTCAGATCACGT", File.replace('R1', 'R2')]) # 16S
	reads_cytB_R2 = subprocess.check_output(["zgrep", "CCCTCAGAATGATATTTGTCCTCA", File.replace('R1', 'R2')]) # cytB
	counter=0
	# identify unique sequences
	haplotypes_16s = makeD(reads_16s)
	# BLAST
	out16=blast(haplotypes_16s[0], haplotypes_16s[1], File, '16S')
	
	haplotypes_16s_R2 = makeD(reads_16s_R2)
	out16R=blast(haplotypes_16s_R2[0], haplotypes_16s_R2[1], File.replace('R1', 'R2'), '16S_r')
	haplotypes_cytB = makeD(reads_cytB)
	outcytB=blast(haplotypes_cytB[0], haplotypes_cytB[1], File, 'cytB')
	haplotypes_cytB_R2 = makeD(reads_cytB_R2)
	outcytB_R=blast(haplotypes_cytB_R2[0], haplotypes_cytB_R2[1], File.replace('R1', 'R2'), 'cytB_r')
	
	
	outlist = [out16,out16R,outcytB,outcytB_R] # a list of all outputs
	for outfile in outlist:
		if outfile is None: # In some cases, there are no BLAST results for a given locus and F/R read.
			continue
		# first I go through all locus and F/R reads and note cases where one top hit was observed
		elif len(outfile)==1:
			species = outfile[0]
			if species not in D:
				D[species]=0
			D[species]+=1
	# Now I go back through the locus and F/R reads and use the previous top hits to be the "tiebreaker" for cases of multiple top hits.
	for outfile in outlist:
		if outfile is None:
			continue
		if len(outfile)>1:
			prev=0
			for prev_hit in D:
				if prev_hit in outfile:
					prev=1
					D[prev_hit]+=1
			# If no 'tiebreaker', then add all new species to list
			if not prev:
				for species in outfile:
					if species not in D:
						D[species]=0
					D[species]+=1
	name=File.split("_")[0]
	print name, "\t".join(["\t".join([SPP,'unique_seq:', str(D[SPP])]) for SPP in D.keys()])
