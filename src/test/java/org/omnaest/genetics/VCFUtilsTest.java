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
package org.omnaest.genetics;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;
import org.omnaest.genetics.domain.VCFData;
import org.omnaest.genetics.domain.VCFRecord;
import org.omnaest.genetics.domain.VCFRecord.SampleFields.GenoType;
import org.omnaest.genetics.domain.VCFRecord.SampleInfo;
import org.omnaest.genetics.translator.domain.CodeAndPosition;
import org.omnaest.genetics.translator.domain.NucleicAcidCode;
import org.omnaest.genetics.translator.domain.NucleicAcidCodeSequence;
import org.omnaest.utils.element.UnaryLeftAndRight;

public class VCFUtilsTest
{

	@Test
	@Ignore
	public void testRead() throws Exception
	{
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example.vcf"))
									.parse();

		vcfData	.getRecords()
				.forEach(record ->
				{
					System.out.println(record);
				});

	}

	@Test
	public void testSampleInfo() throws Exception
	{
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example4.vcf"))
									.parse();

		VCFRecord vcfRecord = vcfData	.getRecords()
										.skip(1)
										.findFirst()
										.get();

		assertEquals("0/1", vcfRecord	.getParsedSampleFields()
										.filterByFieldAsUniqueValue(SampleInfo.GT));
		assertTrue(vcfRecord.getParsedSampleFields()
							.hasGenoType(GenoType.REFERENCE_AND_ALTERNATIVE));
		assertFalse(vcfRecord	.getParsedSampleFields()
								.hasGenoType(GenoType.REFERENCE_BOTH));
		assertFalse(vcfRecord	.getParsedSampleFields()
								.hasGenoType(GenoType.ALTERNATIVE_BOTH));
	}

	@Test
	public void testReadAndStream() throws Exception
	{
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example2.vcf"))
									.parse();

		NucleicAcidCodeSequence referenceSequence = NucleicAcidCodeSequence.valueOf("atCga".toUpperCase());

		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeSequence("X", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atCga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeSequence("1", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atGga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeSequence("2", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeSequence("3", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atCAga".toUpperCase(), mergeSequence);
		}

		{
			AtomicLong position = new AtomicLong(1);
			List<CodeAndPosition<NucleicAcidCode>> mergeSequence = vcfData	.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeCodeAndPositionSequence("3", referenceSequence.stream()
																																			.map(code -> new CodeAndPosition<NucleicAcidCode>(	code,
																																																position.getAndIncrement())))

																			.collect(Collectors.toList());

			int ii = 0;
			assertEquals(1, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("a".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
			assertEquals(2, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("t".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
			assertEquals(3, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("C".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
			assertEquals(4, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("A".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
			assertEquals(5, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("g".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
			assertEquals(6, mergeSequence	.get(ii)
											.getPosition());
			assertEquals("a".toUpperCase(), mergeSequence	.get(ii++)
															.getCode()
															.toString());
		}

	}

	@Test
	public void testApplicator() throws Exception
	{
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example3.vcf"))
									.parse();

		NucleicAcidCodeSequence referenceSequence = NucleicAcidCodeSequence.valueOf("atCga".toUpperCase());

		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.usingPrimaryAllele()
																			.applyToChromosomeSequence("1", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atGga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
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
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example3.vcf"))
									.parse();

		Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> positionToReplacement = vcfData	.applicator()
																							.getPositionToReplacementForChromosome("1");

		assertEquals(1, positionToReplacement.size());
		assertEquals("G", positionToReplacement	.get(3l)
												.get(0)
												.getRight()
												.toString());
		assertEquals("A", positionToReplacement	.get(3l)
												.get(1)
												.getRight()
												.toString());
	}

}
