import vcf
import argparse
from collections import defaultdict
import pandas as pd

# example run:
# python base_script.py -v 1.vcf -v 2.vcf -ac /home/ibab/Applications/semester_2/project/databases/park_samples/00-common_all.vcf -c /home/ibab/Applications/semester_2/project/databases/park_samples/CosmicCodingMuts.vcf -cl /home/ibab/Applications/semester_2/project/databases/park_samples/clinvar.vcf


vcf_variants = defaultdict(set)
vcf_positions = defaultdict(set)
variant_features = defaultdict(dict)
commonMutationsVar = set()
commonMutationPos = set()
clinvarMutations = set()
non_common_var = {}       
non_common_pos = {}
chr_count_common = {}


def parseVCF(vcffile):
	vcfreader = vcf.Reader(open(vcffile))
	samples = vcfreader.samples
	sampleName = None
	features = ['AF', 'AC', 'DP', 'ExcessHet','FS','MQ','QD',
	'GQ', 'DP', 'QUAL']

	if len(samples) == 0:
		print("No sample found, consolidating it as known_mutations")
		sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_known_mutations"
	for record in vcfreader:
		alts = record.ALT
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '')
		if sampleName and "_known_mutations" in sampleName:
			vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
			alts = record.ALT
			keyh = chrm + "_" + str(record.POS)
			variant_features[sampleName][keyh] = {}
			for alt in alts:
				vcf_variants[sampleName].add("_".join([chrm, str(record.POS), record.REF, str(alt)]))
			for field in features:
				if field in ('AC', 'AF'):
					variant_features[sampleName][keyh][field] = record.INFO[field][0]
				elif field in ('DP', 'ExcessHet','FS','MQ','QD'):
					variant_features[sampleName][keyh][field] = record.INFO[field]
				elif field in ('QUAL'):
					variant_features[sampleName][keyh][field] = record.__getattribute__(field)
				else:
					print(field + " not found for " + keyh)
		else:
			for sample in samples:
				sampleName = vcffile.split("/")[-1].replace(".vcf", "") + "_" + "sample"
				try:
				    alts = record.genotype(sample).gt_bases.split("|")
				except:
					continue
				vcf_positions[sampleName].add(chrm + "_" + str(record.POS))
				keyh = chrm + "_" + str(record.POS)
				if len(alts) == 1:
					alts = record.genotype(sample).gt_bases.split("/")
				for alt in alts:
					vcf_variants[sampleName].add("_".join([chrm, str(record.POS), record.REF, str(alt)]))
				variant_features[sampleName][keyh] = {}
				for field in features:
					if field in ('AF', 'AC'):
						variant_features[sampleName][keyh][field] = record.INFO[field][0]
					elif field in ('ExcessHet','FS','MQ','QD') and field in record.INFO:
						variant_features[sampleName][keyh][field] = record.INFO[field]
					elif field in ('GQ', 'DP'):
						variant_features[sampleName][keyh][field] = record.genotype(sample)[field]
					elif field in ('QUAL'):
						variant_features[sampleName][keyh][field] = record.__getattribute__(field)
					#else:
						#print(field + " not found for " + keyh)



##reading dbsnp common snp file:
def parseCommon(common_snp):
	vcfreader = vcf.Reader(open(common_snp))
	#vcfreader._parse_info
	for record in vcfreader: 
		for alt in record.ALT:
			chrm = record.CHROM.replace('chr', '')
			chrm = chrm.replace('Chr', '')
			identifier = "_".join([chrm, str(record.POS), record.REF, str(alt)])
			pos_identifier = "_".join([chrm, str(record.POS)])
			commonMutationsVar.add(identifier)
			commonMutationPos.add(pos_identifier)
			
	return commonMutationsVar, commonMutationPos



##creating variant and position dictionaries of non common mutation
def non_common_mutations(vcf_positions, vcf_variants, commonMutationsVar,commonMutationPos):
	samples = vcf_positions.keys()
	print("SAMPLES being compared:", samples)

	
	for sample1 in samples:
		non_common_var[sample1] = (vcf_variants[sample1] - (commonMutationsVar))
		non_common_pos[sample1] = (vcf_positions[sample1] - (commonMutationPos))
	#for sample1 in samples:
		fout = open(sample1 + "_filteredVariants.txt", 'w')
		for variant in non_common_var[sample1]:
			fout.writelines(variant + "\n")
		fout.close()
		print("No. of variants filtered away : ",sample1, " : ",(len(vcf_variants[sample1]) - len(non_common_var[sample1])))

	return non_common_var,non_common_pos


##creating variant and position dictionaries of cosmic mutations
def parseCosmic(cosmicVCF):
	vcfreader = vcf.Reader(open(cosmicVCF))
	vcfreader._parse_info
	cosmicMutations = set()
	for record in vcfreader: 
		for alt in record.ALT:
			chrm = record.CHROM.replace('chr', '')
			chrm = chrm.replace('Chr', '')
			identifier = "_".join([chrm, str(record.POS), record.REF, str(alt)])
			cosmicMutations.add(identifier)

	return cosmicMutations



##creating variant and position dictionaries of clinvar mutations
def parseClinvar(clinvarVCF):
	vcfreader = vcf.Reader(open(clinvarVCF))
	vcfreader._parse_info
	clinvarMutations = set()
	for record in vcfreader:
		chrm = record.CHROM.replace('chr', '')
		chrm = chrm.replace('Chr', '') 
		for alts in record.ALT:
			alts = str(alts)
			for alt in alts:
				for refs in record.REF:
					identifier = "_".join([chrm, str(record.POS), refs, str(alt)])
					clinvarMutations.add(identifier)
	return clinvarMutations





def common_mutations(non_common_var, non_common_pos, cosmicMutations,clinvarMutations):
	samples = non_common_var.keys()
	print("SAMPLES being compared:", samples)
	
	allCommon = set()
	commonVariants = {}
	cosmicOverlap = {}
	clinvarOverlap = {}
	## Common across all samples:
	for sample in samples:
		if allCommon == set():
			allCommon = non_common_var[sample]
		else:
			allCommon = allCommon.intersection(non_common_var[sample])
	
	fout = open("Common_variants.txt", 'w')
	for variant in allCommon:
		fout.writelines(variant + "\n")
		var = variant.split("_")
		chromo = "chr_"+var[0]
		if var[0].isdigit() or var[0] in "XY":
			if chromo in chr_count_common:
				chr_count_common[chromo] += 1
			else:
				chr_count_common[chromo] = 1
	fout.close()


	print("chromosome-wise counts of common variants: ", chr_count_common)
	print("Common Across all samples:", len(allCommon))

	## Common variants between every 2 samples and with cosmic:
	fileh = open("cosmic_overlap.txt","w")
	fileh2 = open("clinvar_overlap.txt","w")
	
	count1 = 0
	for sample1 in samples:
		count1 += 1
		fileh4 = open(str(count1) + "_samplewise.txt","w")
		commonVariants[sample1] = {}
		count2 = 0
		for sample2 in samples:
			count2 += 1
			if sample1 == sample2:
				continue
			else:
				common_inter = non_common_var[sample1].intersection(non_common_var[sample2])
				commonVariants[sample1][sample2] = len(common_inter)
				if (count1 + 1) == count2:
					for variant in common_inter:
						fileh4.writelines(variant+"\n")
	
					

		co_overlap = (non_common_var[sample1].intersection(cosmicMutations))
		cosmicOverlap[sample1] = co_overlap
		for cos_variants in co_overlap:
			fileh.writelines(cos_variants + "\n")
		print ("COSMIC OVERLAP:", sample1, " : " , len(cosmicOverlap[sample1]))
		cl_overlap = (non_common_var[sample1].intersection(clinvarMutations))
		clinvarOverlap[sample1] = cl_overlap
		for clin_variants in cl_overlap:
			fileh2.writelines(clin_variants + "\n")
		print ("CLINVAR OVERLAP:", sample1, " : ", len(clinvarOverlap[sample1]))
	fileh4.close()
	fileh.close()
	fileh2.close()


	fathmm = {}
	df = pd.read_csv("fathmm_ID.csv",dtype=object)
	cos = df.to_dict("records")
	for key,val in cosmicOverlap.items():
		for value in val:
			values = value.split("_")
			for dicts in cos:
				if values[0] == dicts["chr"] and int(values[1]) >= float(dicts["pos_1"]) and int(values[1]) <= float(dicts["pos_2"]):
					dictval = key+":"+value
					fathmm[dictval] = [dicts["FATHMM_score"],dicts["FATHMM_prediction"],dicts["gene_name"],dicts["Mutation ID"]]



	print("COMMON VARIANTS matix:")
	sample_matrix = pd.DataFrame(commonVariants)
	print (sample_matrix.to_string())
	sample_matrix.to_csv("sample_matrix.csv")
	cos_overlap = pd.DataFrame.from_dict(cosmicOverlap, orient='index')
	cos_overlap.to_csv("cosmic_overlap.csv")
	clin_overlap = pd.DataFrame.from_dict(clinvarOverlap, orient='index')
	clin_overlap.to_csv("clinvar_overlap.csv")
	print("FATHMM SCORES: ", fathmm)
	fathmm_df = pd.DataFrame.from_dict(fathmm, columns= ["Fathmm_score", "Fathmm_prediction","gene_name","Mutation_id"], orient = "index")
	fathmm_df = fathmm_df.T
	#fathmm_df.columns =["Fathmm_score", "Fathmm_prediction"]
	#fathmm_df["Fathmm_score"] = fathmm_df["Fathmm_score"].astype(str)
	fathmm_df.to_csv("fathmm_consolidation.csv")



	fout = open("AllCommon_PCA.txt", 'w')
	line = ["Sample"]
	features = ['AF', 'AC', 'DP', 'ExcessHet','FS','MQ','QD',
	'GQ', 'DP', 'QUAL']
	alreadyseen = set()
	for variant in allCommon: 
		for field in features:
			k = variant + "_" + field
			if k in alreadyseen:
				continue
			else:
				line.append(k)
				alreadyseen.add(k)
	fout.writelines("\t".join(line) + "\n")

	for sample in samples:
		alreadyseen = set()
		line = [sample]
		for variant in allCommon:
			for field in features:
				k = variant + "_" + field
				if k in alreadyseen:
					continue
				else:
					alreadyseen.add(k)
				try:
					if field in variant_features[sample][variant]:
						line.append(str(variant_features[sample][variant][field]))
					else:
						line.append('0')

				except:
					continue
		fout.writelines("\t".join(line) + "\n")
	fout.close()


	return (commonVariants, cosmicOverlap)


def unique_mutations(non_common_var):
	samples = non_common_var.keys()
	vcf_exclusives = {}

	## Unique variants in every sample
	for sample1 in samples:
		chr_count_uniq = {}
		vcf_exclusives[sample1] = non_common_var[sample1]
		for sample2 in samples:
			if sample1 == sample2:
				continue
			else:
				vcf_exclusives[sample1] = vcf_exclusives[sample1] - non_common_var[sample2]

		fout = open(sample1 + "_uniqueVariants.txt", 'w')
		for variant in vcf_exclusives[sample1]:
			fout.writelines(variant + "\n")
			var = variant.split("_")
			chromo = "chr_"+var[0]
			if var[0].isdigit() or var[0] in "XY":
				if chromo in chr_count_uniq:
					chr_count_uniq[chromo] += 1
				else:
					chr_count_uniq[chromo] = 1
		fout.close()

	
		print(sample1 , "chromosome-wise counts of unique variants: " , chr_count_uniq)

	return vcf_exclusives





def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", "--vcfs", action="append", help="list vcfs")
	parser.add_argument("-c", "--cosmic", help="cosmic vcf")
	parser.add_argument("-ac", "--all_common", help="common mutations vcf")
	parser.add_argument("-o", "--output_path", help="output path")
	parser.add_argument("-cl", "--clinvar", help="clinvar vcf")



	
	#parser.add_argument("-o","--output_path", help="output path")

	args = parser.parse_args()

	for eachVcf in args.vcfs:
		print("PARSING", eachVcf)
		parseVCF(eachVcf)


	commonMutations = parseCommon(args.all_common)
	non_common = non_common_mutations(vcf_positions, vcf_variants, commonMutationsVar,commonMutationPos)
	cosmicMutations = parseCosmic(args.cosmic)
	clinvarMutations = parseClinvar(args.clinvar)
	(commonVariants, cosmicOverlap) = common_mutations(non_common_var, non_common_pos, cosmicMutations,clinvarMutations)
	vcf_exclusives = unique_mutations(non_common_var)   
    
	
if __name__ == "__main__":
	main()
