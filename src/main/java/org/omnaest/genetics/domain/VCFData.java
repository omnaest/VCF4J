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
package org.omnaest.genetics.domain;

import java.util.stream.Stream;

import org.omnaest.genetics.fasta.translator.NucleicAcidCode;

public interface VCFData
{

	public static interface GenomeApplicator
	{
		/**
		 * Similar to {@link #usingAllele(int)} with value = 0
		 * 
		 * @return
		 */
		public GenomeApplicator usingPrimaryAllele();

		/**
		 * Similar to {@link #usingAllele(int)} with value = 1
		 * 
		 * @return
		 */
		public GenomeApplicator usingSecondaryAllele();

		/**
		 * If multiple alleles are found for a single position this defines which one is used
		 * 
		 * @param allele
		 *            = 0,1,...
		 * @return
		 */
		public GenomeApplicator usingAllele(int allele);

		/**
		 * Applies the {@link VCFData} to a given {@link Stream} of {@link NucleicAcidCode}s and returns the recombined {@link Stream}.<br>
		 * <br>
		 * This is used to apply
		 * the {@link VCFData} to its reference genome to recalculate the fasta data.
		 * 
		 * @param chromosome
		 * @param sequence
		 * @return
		 */
		public Stream<NucleicAcidCode> applyToChromosome(String chromosome, Stream<NucleicAcidCode> sequence);
	}

	public Stream<VCFRecord> getRecords();

	/**
	 * Returns the {@link GenomeApplicator} instance to apply the {@link VCFRecord}s to its reference genome
	 * 
	 * @return
	 */
	public GenomeApplicator applicator();

}