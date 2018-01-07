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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.omnaest.genetics.components.VCFParserManager;
import org.omnaest.genetics.components.parser.VCFParser;
import org.omnaest.genetics.components.parser.VCFParser_4_1;
import org.omnaest.genetics.domain.VCFData;
import org.omnaest.genetics.domain.VCFRecord;
import org.omnaest.genetics.translator.domain.CodeAndPosition;
import org.omnaest.genetics.translator.domain.NucleicAcidCode;
import org.omnaest.utils.ConsumerUtils;
import org.omnaest.utils.ListUtils;
import org.omnaest.utils.StreamUtils;
import org.omnaest.utils.element.lar.UnaryLeftAndRight;

/**
 * Utils regarding the variant call format<br>
 * <br>
 * Example use:<br>
 * 
 * <pre>
 * VCFData vcfData = VCFUtils.read()
 * 							.from(new File("genome.vcf"))
 * 							.parse();
 * </pre>
 * 
 * @author omnaest
 */
public class VCFUtils
{
	private static VCFParserManager parserManager = new VCFParserManager().register(new VCFParser_4_1());

	public static interface VCFReader
	{
		/**
		 * Similar to {@link #from(File, Charset)} with {@link StandardCharsets#UTF_8}
		 * 
		 * @param file
		 * @return
		 * @throws FileNotFoundException
		 */
		public VCFReader from(File file) throws FileNotFoundException;

		/**
		 * Reads the {@link VCFRecord}s from a {@link File} with the given {@link Charset}
		 * 
		 * @param file
		 * @param charset
		 * @return
		 * @throws FileNotFoundException
		 */
		public VCFReader from(File file, Charset charset) throws FileNotFoundException;

		/**
		 * Reads the {@link VCFRecord}s from an {@link InputStream} using the given {@link Charset}
		 * 
		 * @param inputStream
		 * @param charset
		 * @return
		 */
		public VCFReader from(InputStream inputStream, Charset charset);

		/**
		 * Similar to {@link #from(InputStream, Charset)} with {@link StandardCharsets#UTF_8}
		 * 
		 * @param inputStream
		 * @return
		 */
		public VCFReader from(InputStream inputStream);

		/**
		 * Reads the {@link VCFRecord}s from a given {@link Reader}
		 * 
		 * @param reader
		 * @return
		 */
		public VCFReader from(Reader reader);

		/**
		 * Parses the {@link VCFRecord}s and closes the underlying parser. This operation is not repeatable. This operation does not load the content into
		 * memory an is implemented for large vcf file {@link Stream} processing.
		 * 
		 * @return
		 */
		public Stream<VCFRecord> parseOnce();

		/**
		 * Parses the {@link VCFRecord}s and constructs an in memory {@link VCFData} instance with the complete content
		 * 
		 * @return
		 */
		public VCFData parse();
	}

	public static VCFReader read()
	{
		return new VCFReader()
		{
			private Reader reader;

			@Override
			public VCFReader from(File file) throws FileNotFoundException
			{
				return this.from(file, StandardCharsets.UTF_8);

			}

			@Override
			public VCFReader from(File file, Charset charset) throws FileNotFoundException
			{
				return this.from(new FileInputStream(file), charset);
			}

			@Override
			public VCFReader from(InputStream inputStream, Charset charset)
			{
				return this.from(new InputStreamReader(inputStream, charset));
			}

			@Override
			public VCFReader from(InputStream inputStream)
			{
				return this.from(inputStream, StandardCharsets.UTF_8);
			}

			@Override
			public VCFReader from(Reader reader)
			{
				this.reader = reader;
				return this;
			}

			@Override
			public VCFData parse()
			{
				Map<String, List<VCFRecord>> chromosomeToRecords = this	.parseOnce()
																		.collect(Collectors.groupingBy(record -> StringUtils.upperCase(record.getChromosome())));

				return new VCFData()
				{
					@Override
					public Stream<VCFRecord> getRecords()
					{
						return chromosomeToRecords	.values()
													.stream()
													.flatMap(records -> records.stream());
					}

					@Override
					public GenomeApplicator applicator()
					{
						return new GenomeApplicator()
						{
							@Override
							public AlleleSpecificGenomeApplicator usingSecondaryAllele()
							{
								return this.usingAllele(1);
							}

							@Override
							public AlleleSpecificGenomeApplicator usingPrimaryAllele()
							{
								return this.usingAllele(0);
							}

							@Override
							public Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> getPositionToReplacementForChromosome(String chromosome)
							{
								return this.determinePositionToReplacement(chromosomeToRecords.getOrDefault(StringUtils.upperCase(chromosome),
																											Collections.emptyList()));
							}

							@Override
							public Stream<ChromosomeAndPositionReplacement> getPositionToReplacements()
							{
								return chromosomeToRecords	.keySet()
															.stream()
															.map(chromosome -> new ChromosomeAndPositionReplacement()
															{
																@Override
																public Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> getPositionToReplacement()
																{
																	return getPositionToReplacementForChromosome(chromosome);
																}

																@Override
																public String getChromosome()
																{
																	return chromosome;
																}
															});
							}

							private Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> determinePositionToReplacement(List<VCFRecord> records)
							{
								Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> positionToReplacement = new ConcurrentHashMap<>();

								if (records != null)
								{
									for (VCFRecord vcfRecord : records)
									{
										String reference = vcfRecord.getReference();
										String alternativeAlleles = vcfRecord.getAlternativeAlleles();
										long position = NumberUtils.toLong(vcfRecord.getPosition());

										for (int ii = 0; ii < reference.length() || ii < alternativeAlleles.length(); ii++)
										{
											long currentPosition = position + ii;
											NucleicAcidCode left = ii < reference.length() ? NucleicAcidCode.valueOf(reference.charAt(ii)) : null;
											NucleicAcidCode right = ii < alternativeAlleles.length() ? NucleicAcidCode.valueOf(alternativeAlleles.charAt(ii))
													: null;
											positionToReplacement	.computeIfAbsent(currentPosition, c -> new ArrayList<>())
																	.add(new UnaryLeftAndRight<NucleicAcidCode>(left, right));
										}
									}
								}

								return positionToReplacement;
							}

							@Override
							public AlleleSpecificGenomeApplicator usingAllele(int allele)
							{
								return new AlleleSpecificGenomeApplicator()
								{

									@Override
									public Stream<NucleicAcidCode> applyToChromosomeSequence(String chromosome, Stream<NucleicAcidCode> sequence)
									{
										AtomicLong position = new AtomicLong(1);
										return this	.applyToChromosomeCodeAndPositionSequence(	chromosome,
																								sequence.map(code -> new CodeAndPosition<>(	code,
																																			position.getAndIncrement())))
													.map(cap -> cap.getCode());
									}

									@Override
									public Stream<CodeAndPosition<NucleicAcidCode>> applyToChromosomeCodeAndPositionSequence(	String chromosome,
																																Stream<CodeAndPosition<NucleicAcidCode>> sequence)
									{
										Map<Long, List<UnaryLeftAndRight<NucleicAcidCode>>> positionToReplacement = getPositionToReplacementForChromosome(chromosome);

										AtomicLong position = new AtomicLong(-1);
										return sequence	.peek(ConsumerUtils.consumeOnce(code ->
										{
											if (position.get() < 0)
											{
												position.set(Math.max(1, code.getPosition()));
											}
										}))
														.flatMap(code ->
														{
															//
															Stream<NucleicAcidCode> retval = Stream.of(code.getCode());

															//
															long currentPosition = code.getPosition();
															List<UnaryLeftAndRight<NucleicAcidCode>> replacements = positionToReplacement.get(currentPosition);
															if (replacements != null)
															{
																UnaryLeftAndRight<NucleicAcidCode> replacement = ListUtils.get(replacements, allele);

																if (replacement != null)
																{
																	NucleicAcidCode referenceCode = replacement.getLeft();
																	NucleicAcidCode replacementCode = replacement.getRight();

																	if (referenceCode == null)
																	{
																		retval = Stream.of(replacementCode, code.getCode());
																	}
																	else
																	{
																		//
																		this.assertReferenceCodeMatches(code.getCode(), currentPosition, referenceCode);

																		//
																		if (replacementCode == null)
																		{
																			retval = Stream.empty();
																		}
																		else
																		{
																			retval = Stream.of(replacementCode);
																		}
																	}
																}
															}

															return retval.map(c -> new CodeAndPosition<>(c, position.getAndIncrement()));
														});
									}

									private void assertReferenceCodeMatches(NucleicAcidCode code, long currentPosition, NucleicAcidCode referenceCode)
									{
										if (!code.equals(referenceCode))
										{
											throw new IllegalStateException("Reference code did not match: " + code + "<->" + referenceCode + " at position: "
													+ currentPosition);
										}
									}
								};
							}

						};
					}

				};
			}

			@Override
			public Stream<VCFRecord> parseOnce()
			{
				Stream<String> lines = StreamUtils	.fromReader(this.reader)
													.filter(line -> !StringUtils.isBlank(line));

				VCFParser parser = parserManager.getInstance(lines);

				return parser.getRecords();

			}

		};
	}

	public static VCFParserManager getParserManager()
	{
		return parserManager;
	}

}
