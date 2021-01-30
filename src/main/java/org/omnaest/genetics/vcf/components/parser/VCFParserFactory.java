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
package org.omnaest.genetics.vcf.components.parser;

import java.util.List;
import java.util.stream.Stream;

public interface VCFParserFactory
{
	public interface VCFParserFactoryWithHeader
	{
		public boolean canHandle();

		public VCFParser createInstance(Stream<String> lines);
	}

	public double getVersion();

	public VCFParserFactoryWithHeader withHeaders(List<String> headers);

}