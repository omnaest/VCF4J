/*******************************************************************************
 * Copyright 2021 Danny Kunz
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License.  You may obtain a copy
 * of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
 * License for the specific language governing permissions and limitations under
 * the License.
 ******************************************************************************/
/*

	Copyright 2017 Danny Kunz

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.


*/
package org.omnaest.genomics.vcf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Ignore;
import org.junit.Test;
import org.omnaest.genomics.translator.domain.CodeAndPosition;
import org.omnaest.genomics.translator.domain.NucleicAcidCode;
import org.omnaest.genomics.translator.domain.NucleicAcidCodeSequence;
import org.omnaest.genomics.vcf.VCFUtils;
import org.omnaest.genomics.vcf.domain.VCFData;
import org.omnaest.genomics.vcf.domain.VCFRecord;
import org.omnaest.genomics.vcf.domain.VCFData.Replacements;
import org.omnaest.genomics.vcf.domain.VCFRecord.SampleInfo;
import org.omnaest.genomics.vcf.domain.VCFRecord.SampleFields.Allele;
import org.omnaest.genomics.vcf.domain.VCFRecord.SampleFields.GenoType;

public class VCFUtilsTest
{

    @Test
    @Ignore
    public void testReadRaw() throws Exception
    {
        File file = new File("C:\\Google Drive\\Body odor data\\12 - NA\\raw\\dna\\DX174527_01_var_annotated.vcf");
        Stream<VCFRecord> records = VCFUtils.read()
                                            .from(file)
                                            .parseOnce();
        records.limit(10000)
               .filter(record -> record.getLossOfFunctionPrediction() > 0.1)
               .forEach(record ->
               {
                   System.out.println(record.getLossOfFunctionPrediction());
               });

    }

    @Test
    @Ignore
    public void testReadGenes() throws Exception
    {
        //        File sourceFile = new File("C:\\Z\\data\\6 - my\\raw\\dna\\DX173853_01_var_annotated.vcf");
        //        File targetFile = new File("C:\\\\Z\\\\data\\\\6 - my\\\\raw\\\\dna\\\\DX173853_01_var_annotated_exom.vcf");

        //        File sourceFile = new File("C:\\Z\\data\\4\\raw\\55101705103780_annotated.vcf");
        //        File targetFile = new File("C:\\\\Z\\\\data\\\\4\\\\raw\\\\55101705103780_annotated_SELENBP1.vcf");

        String gene = "FMO3";// "SELENBP1";
        File sourceFile = new File("C:\\Z\\data\\24 - DA\\56001801065129A.snp.ann.vcf");
        File targetFile = new File("C:\\Z\\data\\24 - DA\\56001801065129A.snp.ann_" + gene + ".vcf");

        Stream<VCFRecord> vcfData = VCFUtils.read()
                                            .from(sourceFile)
                                            .parseOnce();

        VCFUtils.write(vcfData.filter(VCFRecord::hasGene)
                              .filter(record ->
                              {
                                  return record.getGene()
                                               .toUpperCase()
                                               .equals(gene);
                              }))
                .into(targetFile);

        //        System.out.println(count);
    }

    @Test
    @Ignore
    public void testRead() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example.vcf"))
                                  .parse();

        vcfData.getRecords()
               .forEach(record ->
               {
                   System.out.println(record);
               });
    }

    @Test
    public void testSampleInfo() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example4.vcf"))
                                  .parse();

        VCFRecord vcfRecord = vcfData.getRecords()
                                     .skip(1)
                                     .findFirst()
                                     .get();

        assertEquals("0/1", vcfRecord.parseSampleFields()
                                     .filterByFieldAsUniqueValue(SampleInfo.GT));
        assertTrue(vcfRecord.parseSampleFields()
                            .hasGenoType(GenoType.REFERENCE_AND_ALTERNATIVE));
        assertFalse(vcfRecord.parseSampleFields()
                             .hasGenoType(GenoType.REFERENCE_BOTH));
        assertFalse(vcfRecord.parseSampleFields()
                             .hasGenoType(GenoType.ALTERNATIVE_BOTH));

        assertEquals(1, vcfRecord.parseSampleFields()
                                 .resolveUniqueAlleleDepth(Allele.REFERENCE));
        assertEquals(3, vcfRecord.parseSampleFields()
                                 .resolveUniqueAlleleDepth(Allele.ALTERNATIVE));
    }

    @Test
    public void testReadAndStream() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example2.vcf"))
                                  .parse();

        NucleicAcidCodeSequence referenceSequence = NucleicAcidCodeSequence.valueOf("atCga".toUpperCase());

        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeSequence("X", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atCga".toUpperCase(), mergeSequence);
        }
        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeSequence("1", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atGga".toUpperCase(), mergeSequence);
        }
        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeSequence("2", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atga".toUpperCase(), mergeSequence);
        }
        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeSequence("3", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atCAga".toUpperCase(), mergeSequence);
        }

        {
            AtomicLong position = new AtomicLong(1);
            List<CodeAndPosition<NucleicAcidCode>> mergeSequence = vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeCodeAndPositionSequence("3", referenceSequence.stream()
                                                                                                                                          .map(code -> new CodeAndPosition<NucleicAcidCode>(code,
                                                                                                                                                                                            position.getAndIncrement())))

                                                                          .collect(Collectors.toList());

            int ii = 0;
            assertEquals(1, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("a".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
            assertEquals(2, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("t".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
            assertEquals(3, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("C".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
            assertEquals(4, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("A".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
            assertEquals(5, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("g".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
            assertEquals(6, mergeSequence.get(ii)
                                         .getPosition());
            assertEquals("a".toUpperCase(), mergeSequence.get(ii++)
                                                         .getCode()
                                                         .toString());
        }

    }

    @Test
    public void testApplicator() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example3.vcf"))
                                  .parse();

        NucleicAcidCodeSequence referenceSequence = NucleicAcidCodeSequence.valueOf("atCga".toUpperCase());

        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingPrimaryAllele()
                                                                          .applyToChromosomeSequence("1", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atGga".toUpperCase(), mergeSequence);
        }
        {
            String mergeSequence = NucleicAcidCodeSequence.valueOf(vcfData.applicator()
                                                                          .usingSecondaryAllele()
                                                                          .applyToChromosomeSequence("1", referenceSequence.stream())
                                                                          .collect(Collectors.toList()))
                                                          .toString();
            assertEquals("atAga".toUpperCase(), mergeSequence);
        }

    }

    @Test
    public void testApplicatorPosition() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example3.vcf"))
                                  .parse();

        Map<Long, Replacements> positionToReplacement = vcfData.applicator()
                                                               .getPositionToReplacementForChromosome("1");

        assertEquals(1, positionToReplacement.size());
        assertEquals("G", positionToReplacement.get(3l)
                                               .getReplacementForAllele(0)
                                               .stream()
                                               .findFirst()
                                               .get()
                                               .getRight()
                                               .toString());
        assertEquals("A", positionToReplacement.get(3l)
                                               .getReplacementForAllele(1)
                                               .stream()
                                               .findFirst()
                                               .get()
                                               .getRight()
                                               .toString());
    }

    @Test
    public void testApplicatorPositionAlleleSpecific() throws Exception
    {
        VCFData vcfData = VCFUtils.read()
                                  .from(this.getClass()
                                            .getResourceAsStream("/example5.vcf"))
                                  .parse();

        Map<Long, Replacements> positionToReplacement = vcfData.applicator()
                                                               .getPositionToReplacementForChromosome("1");

        assertEquals(2, positionToReplacement.size());
        assertEquals("G", positionToReplacement.get(3l)
                                               .getReplacementForAllele(1)
                                               .stream()
                                               .findFirst()
                                               .get()
                                               .getRight()
                                               .toString());
        assertFalse(positionToReplacement.get(3l)
                                         .hasReplacementForAllele(0));
        assertEquals("A", positionToReplacement.get(10l)
                                               .getReplacementForAllele(0)
                                               .iterator()
                                               .next()
                                               .getRight()
                                               .toString());
        assertFalse(positionToReplacement.get(10l)
                                         .hasReplacementForAllele(1));
    }

}
