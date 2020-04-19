package org.omnaest.genetics.domain;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.omnaest.genetics.VCFUtils;

public class VCFRecordTest
{

    @Test
    public void testGetAnnotation() throws Exception
    {
        VCFRecord record = VCFUtils.read()
                                   .from("##fileformat=VCFv4.2\nchr5\t78373251\trs578147669\tG\tGA\t336.0\toff-target\tMQM=60;ANN=GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000255192|protein_coding|1/7|c.34-39dupA||||||INFO_REALIGN_3_PRIME,GA|intron_variant|MODIFIER|DMGDH|ENSG00000132837|transcript|ENST00000520388|processed_transcript|4/4|n.607-21521dupT||||||,GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000518666|protein_coding|2/4|c.-147-39dupA||||||WARNING_TRANSCRIPT_INCOMPLETE&INFO_REALIGN_3_PRIME,GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000521567|protein_coding|1/6|c.34-39dupA||||||INFO_REALIGN_3_PRIME,GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000518758|retained_intron|1/3|n.49-39dupA||||||INFO_REALIGN_3_PRIME,GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000519743|nonsense_mediated_decay|1/6|c.34-39dupA||||||INFO_REALIGN_3_PRIME,GA|intron_variant|MODIFIER|BHMT2|ENSG00000132840|transcript|ENST00000523472|retained_intron|1/1|n.41-39dupA||||||INFO_REALIGN_3_PRIME;GNOMAD_GENOME_AF=0.21182\tGT:DP:AO\t0/1:62:24")
                                   .parseOnce()
                                   .findFirst()
                                   .get();

        assertEquals("BHMT2", record.getAnnotation()
                                    .get("Gene_Name"));
        assertEquals("BHMT2", record.getGene());
    }

    @Test
    public void testLossOfFunctionPrediction() throws Exception
    {
        VCFRecord record = VCFUtils.read()
                                   .from("##fileformat=VCFv4.2\nchr1\t5935162\trs1287637\tA\tT\t2139.0\t.\tMQM=60;OMIM=607215_[NPHP4_(provisional)_Nephronophthisis_4,606966|Senior-Loken_syndrome_4,606996];ANN=T|splice_acceptor_variant&intron_variant|HIGH|NPHP4|ENSG00000131697|transcript|ENST00000378156|protein_coding|20/29|c.2818-2T>A||||||,T|splice_acceptor_variant&intron_variant|HIGH|NPHP4|ENSG00000131697|transcript|ENST00000478423|processed_transcript|16/25|n.2550-2T>A||||||,T|splice_acceptor_variant&intron_variant|HIGH|NPHP4|ENSG00000131697|transcript|ENST00000489180|nonsense_mediated_decay|23/32|c.*629-2T>A||||||,T|splice_acceptor_variant&intron_variant|HIGH|NPHP4|ENSG00000131697|transcript|ENST00000378169|nonsense_mediated_decay|17/26|c.*1719-2T>A||||||,T|splice_acceptor_variant&intron_variant|HIGH|NPHP4|ENSG00000131697|transcript|ENST00000506941|retained_intron|1/1|n.375-2T>A||||||;LOF=(NPHP4|ENSG00000131697|10|0.30);dbNSFP_CADD_phred=18.68;dbNSFP_phyloP100way_vertebrate=0.452000;GNOMAD_AF=0.82766;GNOMAD_GENOME_AF=0.86154;T1000GP_AF=0.843251;EXAC_AC_Hom=35542;EXAC_AF=0.83625;EXAC_AF_AFR=0.89491;EXAC_AF_AMR=0.79441;EXAC_AF_EAS=0.81794;EXAC_AF_NFE=0.83042;EXAC_AF_SAS=0.83054;EXAC_Hom_AFR=3252;EXAC_Hom_NFE=19394;CLINVAR_ACC=RCV000153587.2;CLINVAR_SIG=benign;COSMIC_ID=COSM3751327\tGT:DP:AO\t1/1:69:69")
                                   .parseOnce()
                                   .findFirst()
                                   .get();

        assertEquals(0.3, record.getLossOfFunctionPrediction(), 0.01);
    }

    @Test
    public void testInsertion() throws Exception
    {
        VCFRecord record = VCFUtils.read()
                                   .from("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t55101705103780\n1\t10108\trs62651026\tC\tCC\t9\tnt\tVT=INS;SST=Single;DP=4;DB;DBA=NA;PAF=NA;NP=88;RN=0\tGT:GQ:AB:DP:FC:INS:DEL:AD:QS:LOD:LDF:LDR:PW:PWF:PWR:MMQ:MMQS:IDL:MQ0:MDL:MDR:MDM:CGIANN_VARNAME:CGIANN_1000GAF:CGIANN_ESP6500AF\t0/1:7:C,CC:4:4:3:1:1,3:12,87:0.00,9.18:0.00,4.30:0.00,4.88:0.00,0.85:0.00,0.85:0.00,0.99:16,42:198.0,180.3:0,0:0,0:82,110:68,39:0,1:-,NR_046018.2(DDX11L1) n.-1766_-1765insC:-,NOT_FOUND:-,NOT_FOUND")
                                   .parseOnce()
                                   .findFirst()
                                   .get();

        assertTrue(record.hasInsertion());
        assertFalse(record.hasDeletion());
    }

    @Test
    public void testDeletion() throws Exception
    {
        VCFRecord record = VCFUtils.read()
                                   .from("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t55101705103780\n1\t16362109\trs3085714\tACT\tA\t16\tnt\tVT=DEL;SST=Single;DP=20;DB;DBA=NA;PAF=NA;PC=0.400;RN=3\tGT:GQ:AB:DP:FC:INS:DEL:AD:QS:LOD:LDF:LDR:PW:PWF:PWR:MMQ:MMQS:IDL:MQ0:MDL:MDR:MDM:CGIANN_VARNAME:CGIANN_1000GAF:CGIANN_ESP6500AF\t0/1:176:ACT,A:20:0:0:7:13,7:346,184:0.00,16.12:0.00,0.00:0.00,16.12:0.00,0.97:0.00,0.00:0.00,0.97:60,52:11.1,90.6:0,0:0,0:124,114:27,37:16,10:-,chr1 g.16362110_16362111delCT (intergenic):-,NOT_FOUND:-,NOT_FOUND")
                                   .parseOnce()
                                   .findFirst()
                                   .get();

        assertTrue(record.hasDeletion());
        assertFalse(record.hasInsertion());
    }

}
