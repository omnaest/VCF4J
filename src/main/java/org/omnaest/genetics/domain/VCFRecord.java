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
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

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

	public String getPosition()
	{
		return this.position;
	}

	public String getId()
	{
		return this.id;
	}

	public String getReference()
	{
		return this.reference;
	}

	public String getAlternativeAlleles()
	{
		return this.alternativeAlleles;
	}

	public String getQuality()
	{
		return this.quality;
	}

	public String getFilter()
	{
		return this.filter;
	}

	public String getInfo()
	{
		return this.info;
	}

	public String getFormat()
	{
		return this.format;
	}

	public Map<String, String> getSampleFields()
	{
		return this.sampleFields;
	}

	public Map<String, String> getInfoAsMap()
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

	public String getInfo(AdditionalInfo additionalInfo)
	{
		return this	.getInfoAsMap()
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
