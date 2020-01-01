# compareVCF
compare and analyse any number of given VCF files to get interpretation of variants across them

Process:
    1. A main Python script that takes in sample VCFs, Variant Database files (dbSNP All-Common VCF, COSMIC Coding Mutation VCF, COSMIC mutation details TSV, ClinVar VCF) to generate output files based on the intersections of variants with the variants present in the Database by reducing it into a simpler form.
(dbSNP: variants that are present in more than 5% of population and are not clinically significant is filtered out; COSMIC/ClinVar : to check the presence of reported cancer disease causing variants)
	2.  A trace back Python code to get back the features of extracted variants.
	3.  Annotations and Interpretations from the output files of the code.

Arguments:
    1. Basescript.py 
    • -v: Input VCF(s)
    • -ac: dbSNP 00-common_all.vcf 
    • -c: CosmicCodingMuts.vcf
    • -cl: clinvar.vcf
    • FATHMM.csv file should present in the working directory
    2. Traceback.py
    • -f: Text file for variants to be traced back
    • -v: VCF file(s) with the variant calling information
