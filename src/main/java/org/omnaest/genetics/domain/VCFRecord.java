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

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.omnaest.utils.MapUtils;

public class VCFRecord
{
	private String				chromosome;
	private String				position;
	private String				id;
	private String				reference;
	private String				alternativeAlleles;
	private String				quality;
	private String				filter;
	private String				info;
	private String				format;
	private Map<String, String>	sampleFields;

	public VCFRecord(	String chromosome, String position, String id, String reference, String alternativeAlleles, String quality, String filter, String info,
						String format, Map<String, String> sampleFields)
	{
		super();
		this.chromosome = chromosome;
		this.position = position;
		this.id = id;
		this.reference = reference;
		this.alternativeAlleles = alternativeAlleles;
		this.quality = quality;
		this.filter = filter;
		this.info = info;
		this.format = format;
		this.sampleFields = sampleFields;
	}

	public String getChromosome()
	{
		return this.chromosome;
	}

	/**
	 * Returns the POS column
	 * 
	 * @return
	 */
	public String getPosition()
	{
		return this.position;
	}

	/**
	 * Returns the ID column
	 * 
	 * @return
	 */
	public String getId()
	{
		return this.id;
	}

	/**
	 * Returns the REF column
	 * 
	 * @return
	 */
	public String getReference()
	{
		return this.reference;
	}

	/**
	 * Returns the ALT column
	 * 
	 * @return
	 */
	public String getAlternativeAlleles()
	{
		return this.alternativeAlleles;
	}

	/**
	 * Returns the QUAL column
	 * 
	 * @return
	 */
	public String getQuality()
	{
		return this.quality;
	}

	/**
	 * Returns the FILTER column
	 * 
	 * @return
	 */
	public String getFilter()
	{
		return this.filter;
	}

	/**
	 * Returns the INFO column
	 * 
	 * @see #getInfo(AdditionalInfo)
	 * @see #getParsedInfo()
	 * @return
	 */
	public String getInfo()
	{
		return this.info;
	}

	/**
	 * Returns the format {@link String} like e.g. GT:GQ:AB:DP:FC:INS:DEL:AD:QS
	 * 
	 * @return
	 */
	public String getFormat()
	{
		return this.format;
	}

	/**
	 * Returns the additional sample columns and their {@link String} value. Please consider using {@link #getParsedSampleFields()}.
	 * 
	 * @see #getParsedSampleFields()
	 * @return
	 */
	public Map<String, String> getSampleFields()
	{
		return this.sampleFields;
	}

	public static interface SampleFields
	{
		public static class NonUniqueSampleFieldValueException extends IllegalStateException
		{
			private static final long serialVersionUID = -2326486880712988387L;

			public NonUniqueSampleFieldValueException(String message)
			{
				super(message);
			}
		}

		/**
		 * Returns a {@link Map} of samples and the value of the given {@link SampleInfo} field
		 * 
		 * @param sampleInfo
		 * @return
		 */
		public Map<String, String> filterByField(SampleInfo sampleInfo);

		/**
		 * Returns a {@link List} of all values associated to the given {@link SampleInfo} field
		 * 
		 * @param sampleInfo
		 * @return
		 */
		public List<String> filterByFieldAsValues(SampleInfo sampleInfo);

		/**
		 * Returns the unique value for the given {@link SampleInfo} over all samples. If there is more than one distinct value a
		 * {@link NonUniqueSampleFieldValueException}
		 * is thrown.
		 * 
		 * @param sampleInfo
		 * @return
		 */
		public String filterByFieldAsUniqueValue(SampleInfo sampleInfo);

		/**
		 * Returns the {@link #getSampleFields()} parsed using the given {@link #getFormat()} column
		 * 
		 * @return
		 */
		public Map<String, Map<String, String>> get();

	}

	public SampleFields getParsedSampleFields()
	{
		List<String> formatKeys = org.omnaest.utils.StringUtils	.splitToStream(this.format, ":")
																.collect(Collectors.toList());
		Function<String, Map<String, String>> parsingFunction = value -> MapUtils	.builder()
																					.putAll(formatKeys.stream(),
																							org.omnaest.utils.StringUtils.splitToStream(value, ":"))
																					.build();
		Map<String, Map<String, String>> sampleToSampleFieldToValue = VCFRecord.this.getSampleFields()
																					.entrySet()
																					.stream()
																					.collect(Collectors.toMap(	entry -> entry.getKey(),
																												entry -> parsingFunction.apply(entry.getValue())));

		return new SampleFields()
		{

			@Override
			public Map<String, Map<String, String>> get()
			{
				return sampleToSampleFieldToValue;
			}

			@Override
			public String filterByFieldAsUniqueValue(SampleInfo sampleInfo)
			{
				List<String> sampleFieldValues = this.filterByFieldAsValues(sampleInfo);
				if (sampleFieldValues.size() > 1)
				{
					throw new NonUniqueSampleFieldValueException(sampleInfo + "->" + sampleFieldValues	.stream()
																										.collect(Collectors.joining(",")));
				}
				return sampleFieldValues.stream()
										.findFirst()
										.orElse(null);
			}

			@Override
			public List<String> filterByFieldAsValues(SampleInfo sampleInfo)
			{
				return this	.filterByField(sampleInfo)
							.values()
							.stream()
							.distinct()
							.collect(Collectors.toList());
			}

			@Override
			public Map<String, String> filterByField(SampleInfo sampleInfo)
			{
				return this	.get()
							.entrySet()
							.stream()
							.collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry	.getValue()
																								.get(sampleInfo.toString())));
			}
		};
	}

	/**
	 * Returns the {@link #getInfo()} in a parsed format
	 * 
	 * @see #getInfo(AdditionalInfo)
	 * @return
	 */
	public Map<String, String> getParsedInfo()
	{
		Map<String, String> retmap = new LinkedHashMap<>();

		if (StringUtils.isNotBlank(this.info))
		{
			String[] pairs = StringUtils.split(this.info, ";");
			if (pairs != null)
			{
				for (String pair : pairs)
				{
					String[] keyAndValue = StringUtils.splitPreserveAllTokens(pair, "=");
					if (keyAndValue != null && keyAndValue.length == 2)
					{
						retmap.put(keyAndValue[0], keyAndValue[1]);
					}
				}
			}
		}

		return retmap;
	}

	public enum AdditionalInfo
	{
		Gene, AA, AC, AF, AN, BQ, CIGAR, DB, DP, END, H2, H3, MQ, MQ0, NS, SB, SOMATIC, VALIDATED
	}

	public enum SampleInfo
	{
		/** Order of alleles */
		AB,
		/** Log odds of allele being real vs. due to sequencing error (in order specified by AB */
		LOD,
		IDL,
		/** Number of reads aligned to forward strand */
		FC,
		/** Number of small deletions at this location */
		DEL,
		/** Total Read Depth in Sample */
		DP,
		PW,
		MDL,
		MDM,
		LDR,
		MDR,
		/** Genotype */
		GT,
		AD,
		GQ,
		MQ0,
		MMQS,
		LDF,
		QS,
		PWR,
		MMQ,
		/** Number of small insertions at this location */
		INS,
		PWF
	}

	public String getInfo(AdditionalInfo additionalInfo)
	{
		return this	.getParsedInfo()
					.get(additionalInfo.name());
	}

	@Override
	public String toString()
	{
		return "VCFRecord [chromosome=" + this.chromosome + ", position=" + this.position + ", id=" + this.id + ", reference=" + this.reference
				+ ", alternativeAlleles=" + this.alternativeAlleles + ", quality=" + this.quality + ", filter=" + this.filter + ", info=" + this.info
				+ ", format=" + this.format + ", sampleFields=" + this.sampleFields + "]";
	}

	public boolean hasInfo(AdditionalInfo additionalInfo)
	{
		return this.getInfo(additionalInfo) != null;
	}

}
