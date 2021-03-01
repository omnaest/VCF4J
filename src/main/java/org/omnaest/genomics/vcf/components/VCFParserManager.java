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
package org.omnaest.genomics.vcf.components;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.omnaest.genomics.vcf.components.parser.VCFParser;
import org.omnaest.genomics.vcf.components.parser.VCFParserFactory;
import org.omnaest.utils.StreamUtils;
import org.omnaest.utils.StreamUtils.Drainage;

public class VCFParserManager
{
	private SortedSet<VCFParserFactory> parserFactorys = new TreeSet<>((v1, v2) -> -1 * Double.compare(v1.getVersion(), v2.getVersion()));

	public static class NoParserAvailableException extends RuntimeException
	{
		public NoParserAvailableException(List<String> headers)
		{
			super("No parser found matching headers: " + headers.stream()
																.collect(Collectors.joining("\t")));
		}

		private static final long serialVersionUID = 2605840059087548331L;
	}

	public VCFParserManager register(VCFParserFactory parserFactory)
	{
		this.parserFactorys.add(parserFactory);
		return this;
	}

	/**
	 * @throws NoParserAvailableException
	 *             if no version specific parser is available
	 * @param lines
	 * @return
	 */
	public VCFParser getInstance(Stream<String> lines)
	{
		Drainage<String> drainage = StreamUtils.drain(lines, line -> !StringUtils.startsWith(line, "#"));

		List<String> headers = drainage	.getPrefetch()
										.collect(Collectors.toList());
		return this.parserFactorys	.stream()
									.map(parser -> parser.withHeaders(headers))
									.filter(parser -> parser.canHandle())
									.findFirst()
									.orElseThrow(() -> new NoParserAvailableException(headers))
									.createInstance(drainage.getStreamIncludingPrefetch());
	}

}
