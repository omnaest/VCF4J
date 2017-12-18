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

import java.util.stream.Collectors;

import org.junit.Ignore;
import org.junit.Test;
import org.omnaest.genetics.domain.VCFData;
import org.omnaest.genetics.fasta.domain.NucleicAcidCodeSequence;

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
	public void testReadAndStream() throws Exception
	{
		VCFData vcfData = VCFUtils	.read()
									.from(this	.getClass()
												.getResourceAsStream("/example2.vcf"))
									.parse();

		NucleicAcidCodeSequence referenceSequence = NucleicAcidCodeSequence.valueOf("atCga".toUpperCase());

		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.applyToChromosome("X", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atCga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.applyToChromosome("1", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atGga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.applyToChromosome("2", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atga".toUpperCase(), mergeSequence);
		}
		{
			String mergeSequence = NucleicAcidCodeSequence	.valueOf(vcfData.applicator()
																			.applyToChromosome("3", referenceSequence.stream())
																			.collect(Collectors.toList()))
															.toString();
			assertEquals("atCAga".toUpperCase(), mergeSequence);
		}

	}

}
