# VCF4J
Java utils for parsing the variant call format

# Example

 VCFData vcfData = VCFUtils.read()
						   .from(new File("genome.vcf"))
  						   .parse();
