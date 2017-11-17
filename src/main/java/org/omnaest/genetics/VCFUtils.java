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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.omnaest.genetics.components.VCFParserManager;
import org.omnaest.genetics.components.VCFParserManager.NoParserAvailableException;
import org.omnaest.genetics.components.parser.VCFParser;
import org.omnaest.genetics.components.parser.VCFParser_4_1;
import org.omnaest.genetics.domain.VCFData;
import org.omnaest.genetics.domain.VCFRecord;
import org.omnaest.utils.StreamUtils;

public class VCFUtils
{
	private static VCFParserManager parserManager = new VCFParserManager().register(new VCFParser_4_1());

	/**
	 * Similar to {@link #read(File, Charset)} using {@link StandardCharsets#UTF_8}
	 * 
	 * @param file
	 * @return
	 * @throws FileNotFoundException
	 * @throws NoParserAvailableException
	 */
	public static VCFData read(File file) throws FileNotFoundException, NoParserAvailableException
	{
		return read(file, StandardCharsets.UTF_8);
	}

	/**
	 * Similar to {@link #read(InputStream, Charset)} using the given {@link File}
	 * 
	 * @param file
	 * @param charset
	 * @return
	 * @throws FileNotFoundException
	 * @throws NoParserAvailableException
	 */
	public static VCFData read(File file, Charset charset) throws FileNotFoundException, NoParserAvailableException
	{
		return read(new FileInputStream(file), charset);
	}

	/**
	 * Similar to {@link #read(InputStream, Charset)} with {@link StandardCharsets#UTF_8}
	 * 
	 * @param inputStream
	 * @return
	 * @throws NoParserAvailableException
	 */
	public static VCFData read(InputStream inputStream) throws NoParserAvailableException
	{
		return read(inputStream, StandardCharsets.UTF_8);
	}

	/**
	 * Similar to {@link #read(Reader)} using the given {@link InputStream} and {@link Charset}
	 * 
	 * @param inputStream
	 * @param charset
	 * @return
	 * @throws NoParserAvailableException
	 */
	public static VCFData read(InputStream inputStream, Charset charset) throws NoParserAvailableException
	{
		return read(new InputStreamReader(inputStream, charset));
	}

	/**
	 * Reads the content of the given {@link Reader} to {@link VCFData} instance
	 * 
	 * @param reader
	 * @return
	 * @throws NoParserAvailableException
	 */
	public static VCFData read(Reader reader) throws NoParserAvailableException
	{
		Stream<String> lines = StreamUtils	.fromReader(reader)
											.filter(line -> !StringUtils.isBlank(line));

		VCFParser parser = parserManager.getInstance(lines);
		return new VCFData()
		{
			@Override
			public Stream<VCFRecord> getRecords()
			{
				return parser.getRecords();
			}
		};
	}

	public static VCFParserManager getParserManager()
	{
		return parserManager;
	}

}
