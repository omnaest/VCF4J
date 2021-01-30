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
package org.omnaest.genetics.vcf;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.io.Writer;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.output.FileWriterWithEncoding;
import org.apache.commons.lang3.StringUtils;
import org.omnaest.genetics.vcf.components.GenomeApplicatorImpl;
import org.omnaest.genetics.vcf.components.VCFParserManager;
import org.omnaest.genetics.vcf.components.parser.VCFParser;
import org.omnaest.genetics.vcf.components.parser.VCFParser_4_1;
import org.omnaest.genetics.vcf.domain.VCFData;
import org.omnaest.genetics.vcf.domain.VCFRecord;
import org.omnaest.utils.IterableUtils;
import org.omnaest.utils.ListUtils;
import org.omnaest.utils.MatcherUtils;
import org.omnaest.utils.MatcherUtils.Match;
import org.omnaest.utils.PatternUtils;
import org.omnaest.utils.StreamUtils;
import org.omnaest.utils.element.bi.BiElement;

/**
 * Utils regarding the variant call format<br>
 * <br>
 * Example use:<br>
 * 
 * <pre>
 * VCFData vcfData = VCFUtils.read()
 *                           .from(new File("genome.vcf"))
 *                           .parse();
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
         * Similar to {@link #from(InputStream)}
         * 
         * @param data
         * @return
         */
        public VCFReader from(byte[] data);

        /**
         * Reads the {@link VCFRecord}s from a given {@link Reader}
         * 
         * @param reader
         * @return
         */
        public VCFReader from(Reader reader);

        /**
         * Reads the {@link VCFRecord}s from a given {@link String}
         * 
         * @param vcfContent
         * @return
         */
        public VCFReader from(String vcfContent);

        /**
         * Parses the {@link VCFRecord}s and closes the underlying parser. This operation is not repeatable. This operation does not load the content into
         * memory and is implemented for large vcf file {@link Stream} processing.
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
            public VCFReader from(byte[] data)
            {
                return this.from(new ByteArrayInputStream(data));
            }

            @Override
            public VCFReader from(Reader reader)
            {
                this.reader = IOUtils.toBufferedReader(reader, 32 * 1024 * 1024);
                return this;
            }

            @Override
            public VCFReader from(String vcfContent)
            {
                return this.from(new StringReader(vcfContent));
            }

            @Override
            public VCFData parse()
            {
                BiElement<Stream<VCFRecord>, Map<String, List<String>>> recordsWithComments = this.parseOnceWithComments();
                Map<String, List<VCFRecord>> chromosomeToRecords = recordsWithComments.getFirst()
                                                                                      .collect(Collectors.groupingBy(record -> StringUtils.replaceAll(StringUtils.upperCase(record.getChromosome()),
                                                                                                                                                      "CHR",
                                                                                                                                                      "")));
                Map<String, List<String>> comments = recordsWithComments.getSecond();

                return new VCFData()
                {
                    @Override
                    public Stream<VCFRecord> getRecords()
                    {
                        return chromosomeToRecords.values()
                                                  .stream()
                                                  .flatMap(records -> records.stream());
                    }

                    @Override
                    public GenomeApplicator applicator()
                    {
                        return new GenomeApplicatorImpl(chromosomeToRecords);
                    }

                    @Override
                    public VCFMetaInfo getMetaInfo()
                    {
                        return new VCFMetaInfo()
                        {
                            @Override
                            public String getReference()
                            {
                                return ListUtils.first(comments.get("reference"));
                            }

                            @Override
                            public String getParsedHumanReferenceGenome()
                            {
                                String reference = this.getReference();
                                Optional<Stream<Match>> matches = PatternUtils.matcher()
                                                                              .of(Pattern.compile("hg[0-9]+|GRCH[0-9]+", Pattern.CASE_INSENSITIVE))
                                                                              .findIn(reference);
                                if (!matches.isPresent())
                                {
                                    return null;
                                }
                                else
                                {
                                    return matches.get()
                                                  .findFirst()
                                                  .get()
                                                  .getMatchRegion();
                                }
                            }

                            @Override
                            public String getFileFormat()
                            {
                                return ListUtils.first(comments.get("fileformat"));
                            }

                            @Override
                            public String getFileDate()
                            {
                                return ListUtils.first(comments.get("fileDate"));
                            }

                            @Override
                            public SampleInfos getSampleInfos()
                            {
                                Map<String, Map<String, String>> retmap = new LinkedHashMap<>();

                                String sampleStr = ListUtils.first(comments.get("SAMPLE"));

                                MatcherUtils.matcher()
                                            .of(Pattern.compile("\\<([^\\>]*)\\>"))
                                            .findIn(sampleStr)
                                            .ifPresent(matches ->
                                            {
                                                matches.forEach(match ->
                                                {
                                                    String singleSampleStr = match.getSubGroupsAsStream()
                                                                                  .findFirst()
                                                                                  .orElse(null);
                                                    Map<String, String> sampleMap = new LinkedHashMap<>();
                                                    org.omnaest.utils.StringUtils.splitToStream(singleSampleStr, ",")
                                                                                 .forEach(keyAndValue ->
                                                                                 {
                                                                                     MatcherUtils.matcher()
                                                                                                 .of(Pattern.compile("([^\\=]+)\\=(.*)"))
                                                                                                 .matchAgainst(keyAndValue)
                                                                                                 .map(keyAndValueMatch -> keyAndValueMatch.getGroups())
                                                                                                 .ifPresent(keyAndValueMatchGroups ->
                                                                                                 {
                                                                                                     String key = keyAndValueMatchGroups.get(1);
                                                                                                     String value = keyAndValueMatchGroups.get(2);

                                                                                                     sampleMap.put(key, value);
                                                                                                 });
                                                                                 });

                                                    String id = sampleMap.get("ID");
                                                    retmap.put(id, sampleMap);
                                                });
                                            });

                                return new SampleInfos()
                                {
                                    @Override
                                    public Map<String, String> getSampleInfo(String id)
                                    {
                                        return retmap.get(id);
                                    }

                                    @Override
                                    public Set<String> getIds()
                                    {
                                        return retmap.keySet();
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
                return this.parseOnceWithComments()
                           .getFirst();

            }

            public BiElement<Stream<VCFRecord>, Map<String, List<String>>> parseOnceWithComments()
            {
                Stream<String> lines = StreamUtils.fromReaderAsLines(this.reader)
                                                  .filter(line -> !StringUtils.isBlank(line));

                VCFParser parser = parserManager.getInstance(lines);

                return BiElement.of(parser.getRecords(), parser.getComments());

            }

        };
    }

    public static VCFParserManager getParserManager()
    {
        return parserManager;
    }

    public static interface VCFWriter
    {

        void into(Writer writer) throws IOException;

        void into(File file) throws IOException;

        void into(File file, Charset encoding) throws IOException;

    }

    public static VCFWriter write(Stream<VCFRecord> vcfData)
    {
        return new VCFWriter()
        {
            @Override
            public void into(Writer writer) throws IOException
            {
                //
                writer.write("##fileformat=VCFv4.3\n");
                writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t000000001\n");

                //
                for (VCFRecord record : IterableUtils.from(vcfData.iterator()))
                {
                    List<String> list = new ArrayList<>();

                    list.add(record.getChromosome());
                    list.add(record.getPosition());
                    list.add(record.getId());
                    list.add(record.getReference());
                    list.add(record.getAlternativeAlleles());
                    list.add(record.getQuality());
                    list.add(record.getFilter());
                    list.add(record.getInfo());
                    list.add(record.getFormat());
                    list.addAll(record.getSampleFields()
                                      .values());

                    writer.write(list.stream()
                                     .collect(Collectors.joining("\t")));
                    writer.write("\n");
                }

                //
                writer.flush();
                writer.close();
            }

            @Override
            public void into(File file, Charset encoding) throws IOException
            {
                FileUtils.forceMkdirParent(file);
                this.into(new FileWriterWithEncoding(file, encoding));
            }

            @Override
            public void into(File file) throws IOException
            {
                this.into(file, StandardCharsets.UTF_8);
            }
        };
    }

}
