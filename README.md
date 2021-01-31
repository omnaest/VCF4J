# VCF4J

Java utility for parsing the variant call format file content

See [Variant Call Format Wikipedia](https://en.wikipedia.org/wiki/Variant_Call_Format)

# Example
## Parsing

     VCFData vcfData = VCFUtils.read()
                               .from(new File("genome.vcf"))
                               .parse();

## Parsing and writing

    VCFUtils.write(VCFUtils.read()
                           .fromFile("input.vcf")
                           .parseOnce()
                           .filter(record -> record.getInfoValue(AdditionalInfo.CLIN_risk_factor)
                                                   .isPresent()))
            .intoFile("output.vcf");
            
# Maven Snapshots

    <dependency>
      <groupId>org.omnaest.genomics</groupId>
      <artifactId>VCF4J</artifactId>
      <version>0.0.1-SNAPSHOT</version>
    </dependency>
    
    <repositories>
        <repository>
            <id>ossrh</id>
            <url>https://oss.sonatype.org/content/repositories/snapshots</url>
            <snapshots>
                <enabled>true</enabled>
            </snapshots>
        </repository>
    </repositories>