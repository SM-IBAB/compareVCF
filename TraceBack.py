import vcf
import argparse

def get_back(each_file,vcf_file):
	fh = open(each_file, "r")
	chrs = []
	for lines in fh:
		variants = lines.split("_")
		if variants[1].isnumeric() == False:
			if len(variants) <= 2:
				if variants[2].isnumeric() == False:
					a = "_".join([variants[0],variants[1],variants[2]])
					b = variants[3]
				else:
					a = "_".join([variants[0],variants[1]])
					b = variants[2]
		else:
			a = variants[0]
			b = variants[1]
		chrs.append([a,b])
	vcfFile = vcf.Reader(open(vcf_file))
	newfile = str(each_file) + str(vcf_file) + "_traceBack.vcf"
	with open(newfile, "w") as outfile:
		writer = vcf.Writer(outfile, vcfFile)
		written_variants = 0
		for record in vcfFile:
			for a,b in chrs:
				try:
					if "chr" in (record.CHROM).lower():
						a = ("chr" + a)
					else:
						a = a
					if a in (record.CHROM).lower() and record.POS == int(b):
						writer.write_record(record)
						written_variants += 1
				except:
					continue
		print(written_variants)
					
				
	

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--filtered_variants", action="append", help="filtered_variant_files(output of noncommon function)")
	parser.add_argument("-v", "--vcfs", action="append", help="list vcfs")
	parser.add_argument("-o","--output_path", help="output path")

	args = parser.parse_args()
	count=0
	for eachVcf,each_filt_file in zip(args.vcfs,args.filtered_variants):
		print("PARSING", each_filt_file)
		get_back(each_filt_file,eachVcf)
if __name__ == "__main__":
    main()
