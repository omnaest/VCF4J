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
package org.omnaest.genomics.vcf.components.parser;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.omnaest.genomics.vcf.domain.VCFRecord;
import org.omnaest.utils.MatcherUtils;
import org.omnaest.utils.MatcherUtils.Match;
import org.omnaest.utils.PatternUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class VCFParser_4_1 implements VCFParserFactory
{
    private static final Logger LOG = LoggerFactory.getLogger(VCFParser_4_1.class);

    protected static class CommentFilter implements Predicate<String>
    {
        private Map<String, List<String>> commentMap = new ConcurrentHashMap<>();

        @Override
        public boolean test(String line)
        {
            boolean isCommentLine = StringUtils.startsWith(line, "#");
            if (isCommentLine)
            {
                Pattern pattern = Pattern.compile("\\#\\#([a-zA-Z0-9]+)\\=(.*)");
                MatcherUtils.matcher()
                            .of(pattern)
                            .matchAgainst(line)
                            .map(Match::getGroups)
                            .ifPresent(match ->
                            {
                                this.commentMap.computeIfAbsent(match.get(1), k -> new ArrayList<>())
                                               .add(match.get(2));
                            });
            }
            return !isCommentLine;
        }

        public Map<String, List<String>> getCommentMap()
        {
            return this.commentMap;
        }

    }

    public VCFParser_4_1()
    {
        super();
    }

    @Override
    public double getVersion()
    {
        return 4.1;
    }

    @Override
    public VCFParserFactoryWithHeader withHeaders(List<String> headers)
    {
        String version = this.determineVersion(headers);
        Map<Integer, String> columnIndexToField = this.determineColumns(headers);
        return new VCFParserFactoryWithHeader()
        {
            @Override
            public boolean canHandle()
            {
                return version != null;
            }

            @Override
            public VCFParser createInstance(Stream<String> lines)
            {
                return new VCFParser()
                {
                    private CommentFilter commentFilter = new CommentFilter();

                    @Override
                    public Stream<VCFRecord> getRecords()
                    {
                        return lines.filter(this.commentFilter)
                                    .map(this::mapToRecord);
                    }

                    @Override
                    public Map<String, List<String>> getComments()
                    {
                        return this.commentFilter.getCommentMap();
                    }

                    private VCFRecord mapToRecord(String line)
                    {
                        Map<String, String> retmap = new LinkedHashMap<>();
                        String[] tokens = StringUtils.splitPreserveAllTokens(line, "\t");
                        for (int ii = 0; ii < tokens.length; ii++)
                        {
                            String field = columnIndexToField.get(ii);
                            String value = tokens[ii];
                            if (StringUtils.isBlank(field))
                            {
                                LOG.warn("Unmapped field value: " + value + "(" + line + ")");
                            }
                            else
                            {
                                retmap.put(field, value);
                            }
                        }

                        //
                        String chromosome = retmap.remove("CHROM");
                        String position = retmap.remove("POS");
                        String id = retmap.remove("ID");
                        String reference = retmap.remove("REF");
                        String alternativeAlleles = retmap.remove("ALT");
                        String quality = retmap.remove("QUAL");
                        String filter = retmap.remove("FILTER");
                        String info = retmap.remove("INFO");
                        String format = retmap.remove("FORMAT");
                        Map<String, String> sampleFields = retmap;
                        return new VCFRecord(chromosome, position, id, reference, alternativeAlleles, quality, filter, info, format, sampleFields);
                    }
                };
            }
        };
    }

    private Map<Integer, String> determineColumns(List<String> headers)
    {
        Map<Integer, String> retmap = new LinkedHashMap<>();

        String columnLine = headers.stream()
                                   .filter(line -> StringUtils.startsWith(line, "#") && !StringUtils.startsWith(line, "##"))
                                   .findFirst()
                                   .orElse(null);
        if (columnLine == null)
        {
            columnLine = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        }

        String[] tokens = StringUtils.splitPreserveAllTokens(StringUtils.removeStart(columnLine, "#")
                                                                        .trim(),
                                                             "\t");
        for (int ii = 0; ii < tokens.length; ii++)
        {
            retmap.put(ii, tokens[ii]);
        }

        return retmap;
    }

    private String determineVersion(List<String> headers)
    {
        Pattern pattern = Pattern.compile("[\\#]+fileformat\\=VCFv(4\\..)", Pattern.CASE_INSENSITIVE);
        String version = headers.stream()
                                .map(line -> PatternUtils.matchToGroups(pattern, line))
                                .filter(map -> !map.isEmpty())
                                .map(map -> map.get(1))
                                .findFirst()
                                .orElse(null);
        return version;
    }

}
