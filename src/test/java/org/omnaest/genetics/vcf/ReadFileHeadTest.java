package org.omnaest.genetics.vcf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.junit.Ignore;
import org.junit.Test;
import org.omnaest.genetics.vcf.VCFUtils;
import org.omnaest.genetics.vcf.domain.VCFRecord;
import org.omnaest.genetics.vcf.domain.VCFRecord.AdditionalInfo;
import org.omnaest.utils.FileUtils;

public class ReadFileHeadTest
{
    @Test
    @Ignore
    public void test() throws IOException
    {
        //        String fileName = "C:\\Z\\data\\6 - my\\raw\\dna\\DX173853_01_var_annotated.vcf";
        //        String fileName = "C:\\Z\\data\\25 - GD\\dna\\15001711232461A.snp.dbsnp.ann.vcf";
        //"C:\\Z\\data\\24 - DA\\56001801065129A.snp.ann.vcf";
        // "C:\\Z\\data\\24 - DA\\56001801065129A.snp.vcf";
        String fileName = "C:\\Z\\databases\\ensembl\\current\\variation\\homo_sapiens_incl_consequences-chr13.vcf";
        File file = new File(fileName);

        BufferedWriter writer = new BufferedWriter(new FileWriter(new File(fileName + ".head")));

        FileUtils.read()
                 .from(file)
                 .getAsLinesStream()
                 .limit(1000000)
                 .forEach(line ->
                 {
                     try
                     {
                         writer.append(line);
                         writer.newLine();
                     }
                     catch (IOException e)
                     {
                         throw new IllegalStateException(e);
                     }
                 });

        writer.close();
    }

    @Test
    @Ignore
    public void testVCFPreview() throws IOException
    {
        //        String fileName = "C:\\Z\\data\\6 - my\\raw\\dna\\DX173853_01_var_annotated.vcf";
        //        String fileName = "C:\\Z\\data\\25 - GD\\dna\\15001711232461A.snp.dbsnp.ann.vcf";
        //"C:\\Z\\data\\24 - DA\\56001801065129A.snp.ann.vcf";
        // "C:\\Z\\data\\24 - DA\\56001801065129A.snp.vcf";
        String fileName = "C:\\Z\\databases\\ensembl\\current\\variation\\homo_sapiens_incl_consequences-chr13.vcf";
        File file = new File(fileName);

        Stream<VCFRecord> records = VCFUtils.read()
                                            .from(file)
                                            .parseOnce();
        List<VCFRecord> vcfData = records.filter(record -> StringUtils.equalsAny(record.getInfoTokensValue(AdditionalInfo.VE, 0)
                                                                                       .orElse(null),
                                                                                 "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
                                                                                 "intron_variant", "regulatory_region_variant", "upstream_gene_variant",
                                                                                 "downstream_gene_variant"))
                                         //                                         .limit(10)
                                         .peek(record ->
                                         {
                                             System.out.println(record);
                                             String ve = record.getInfoTokensValue(AdditionalInfo.VE, 0)
                                                               .orElse(null);
                                             System.out.println(ve);
                                         })
                                         .collect(Collectors.toList());

        VCFUtils.write(vcfData.stream())
                .into(new File(file.getAbsolutePath() + "_filtered_relevant.vcf"));
    }

    @Test
    @Ignore
    public void testFilterRecordsWithoutDBSNPIds() throws IOException
    {
        File file = new File("C:\\Google Drive\\Body odor data\\12 - NA\\raw\\dna\\DX174527_01_var_annotated.vcf");
        Stream<VCFRecord> records = VCFUtils.read()
                                            .from(file)
                                            .parseOnce();
        List<VCFRecord> vcfData = records
                                         //        .limit(10000)
                                         .filter(record -> StringUtils.equals(".", record.getId()))
                                         //               .map(record -> record.getGene())
                                         //                                         .distinct()
                                         //                                         .sorted()
                                         .peek(record ->
                                         {
                                             System.out.println(record);
                                         })
                                         .collect(Collectors.toList());

        VCFUtils.write(vcfData.stream())
                .into(new File(file.getAbsolutePath() + "_noRSID.vcf"));
    }
}
