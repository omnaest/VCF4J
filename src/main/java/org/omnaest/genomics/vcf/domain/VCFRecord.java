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
package org.omnaest.genomics.vcf.domain;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.lang3.StringUtils;
import org.omnaest.utils.CollectorUtils;
import org.omnaest.utils.ListUtils;
import org.omnaest.utils.MapUtils;
import org.omnaest.utils.StreamUtils;

/**
 * Representation of a single record within a VCF file
 * 
 * @see #parseInfo()
 * @see #parseSampleFields()
 * @author omnaest
 */
public class VCFRecord
{
    public static final String SEMICOLON = ";";
    public static final String DOT       = ".";

    private String              chromosome;
    private String              position;
    private String              id;
    private String              reference;
    private String              alternativeAlleles;
    private String              quality;
    private String              filter;
    private String              info;
    private String              format;
    private Map<String, String> sampleFields;

    public VCFRecord(String chromosome, String position, String id, String reference, String alternativeAlleles, String quality, String filter, String info,
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
     * Returns a instance of a {@link Set} which returns all by {@value #SEMICOLON} separated values of the {@link #getId()} field. This ignores also
     * {@value #DOT} values.
     * 
     * @return
     */
    public Set<String> getIds()
    {
        return Optional.ofNullable(this.id)
                       .map(id -> StringUtils.splitByWholeSeparator(id, SEMICOLON))
                       .map(Arrays::asList)
                       .map(List::stream)
                       .orElse(Stream.empty())
                       .map(StringUtils::trim)
                       .filter(StringUtils::isNotBlank)
                       .filter(id -> !StringUtils.equals(DOT, id))
                       .collect(Collectors.toSet());
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
     * @see #parseInfo()
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
     * Returns the additional sample columns and their {@link String} value. Please consider using {@link #parseSampleFields()}.
     * 
     * @see #parseSampleFields()
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

        public enum GenoType
        {
            /** 0/0 or 0 */
            REFERENCE_BOTH("0/0", "0"),
            /** 1/0 or 0/1 */
            REFERENCE_AND_ALTERNATIVE("1/0", "0/1"),
            /** 1/1 or 1 */
            ALTERNATIVE_BOTH("1", "1/1"),
            /**
             * 1/0
             */
            FIRST_ALLELE_ALTERNATIVE("1/0"),
            /**
             * 0/1
             */
            SECOND_ALLELE_ALTERNATIVE("0/1");

            private String[] matchingCodes;

            private GenoType(String... matchingCodes)
            {
                this.matchingCodes = matchingCodes;
            }

            public boolean matches(String code)
            {
                return Arrays.asList(matchingCodes)
                             .stream()
                             .anyMatch(icode -> StringUtils.equalsIgnoreCase(icode, code));
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

        /**
         * Returns true if any {@link #filterByFieldAsValues(SampleInfo)} contains the given {@link GenoType}
         * 
         * @param genoType
         * @return
         */
        public boolean hasGenoType(GenoType genoType);

        public enum Allele
        {
            REFERENCE, ALTERNATIVE
        }

        /**
         * Returns the {@link Allele} depths (AD) if a unique sample is available
         * 
         * @return
         */
        public int resolveUniqueAlleleDepth(Allele allele);

        /**
         * Returns the coverage depth (DP)
         * 
         * @return
         */
        public int resolveUniqueCoverageDepth();

    }

    /**
     * Parses the {@link #getSampleFields()}
     * 
     * @return
     */
    public SampleFields parseSampleFields()
    {
        List<String> formatKeys = org.omnaest.utils.StringUtils.splitToStream(this.format, ":")
                                                               .collect(Collectors.toList());
        Function<String, Map<String, String>> parsingFunction = value -> MapUtils.builder()
                                                                                 .putAll(formatKeys.stream(),
                                                                                         org.omnaest.utils.StringUtils.splitToStream(value, ":"))
                                                                                 .build();
        Map<String, Map<String, String>> sampleToSampleFieldToValue = VCFRecord.this.getSampleFields()
                                                                                    .entrySet()
                                                                                    .stream()
                                                                                    .collect(Collectors.toMap(entry -> entry.getKey(),
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
                    throw new NonUniqueSampleFieldValueException(sampleInfo + "->" + sampleFieldValues.stream()
                                                                                                      .collect(Collectors.joining(",")));
                }
                return sampleFieldValues.stream()
                                        .findFirst()
                                        .orElse(null);
            }

            @Override
            public List<String> filterByFieldAsValues(SampleInfo sampleInfo)
            {
                return this.filterByField(sampleInfo)
                           .values()
                           .stream()
                           .distinct()
                           .collect(Collectors.toList());
            }

            @Override
            public boolean hasGenoType(GenoType genoType)
            {
                return this.filterByFieldAsValues(SampleInfo.GT)
                           .stream()
                           .anyMatch(code -> genoType.matches(code));
            }

            @Override
            public Map<String, String> filterByField(SampleInfo sampleInfo)
            {
                return this.get()
                           .entrySet()
                           .stream()
                           .filter(entry -> entry != null && entry.getKey() != null && entry.getValue() != null && entry.getValue()
                                                                                                                        .get(sampleInfo.toString()) != null)
                           .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()
                                                                                            .get(sampleInfo.toString())));
            }

            @Override
            public int resolveUniqueAlleleDepth(Allele allele)
            {
                String alleleCode = Allele.REFERENCE.equals(allele) ? VCFRecord.this.getReference() : VCFRecord.this.getAlternativeAlleles();
                return this.getAlleleToUniqueAlleleDepths()
                           .getOrDefault(alleleCode, 0);
            }

            @Override
            public int resolveUniqueCoverageDepth()
            {
                return NumberUtils.toInt(this.filterByFieldAsUniqueValue(SampleInfo.DP));
            }

            private Map<String, Integer> getAlleleToUniqueAlleleDepths()
            {
                List<Integer> values = org.omnaest.utils.StringUtils.splitToStream(this.filterByFieldAsUniqueValue(SampleInfo.AD), ",")
                                                                    .map(value -> NumberUtils.toInt(value))
                                                                    .collect(Collectors.toList());
                List<String> keys = org.omnaest.utils.StringUtils.splitToStream(this.filterByFieldAsUniqueValue(SampleInfo.AB), ",")
                                                                 .collect(Collectors.toList());

                return MapUtils.builder()
                               .putAll(keys, values)
                               .build();
            }

        };
    }

    /**
     * Returns the {@link #getInfo()} in a parsed format
     * 
     * @see #getInfo(AdditionalInfo)
     * @return
     */
    public Map<String, String> parseInfo()
    {
        Map<String, String> retmap = new LinkedHashMap<>();

        if (StringUtils.isNotBlank(this.info))
        {
            String[] pairs = StringUtils.split(this.info, SEMICOLON);
            if (pairs != null)
            {
                for (String pair : pairs)
                {
                    String[] keyAndValue = StringUtils.splitPreserveAllTokens(pair, "=");
                    if (keyAndValue != null && keyAndValue.length >= 1)
                    {
                        retmap.put(keyAndValue[0], keyAndValue.length >= 2 ? keyAndValue[1] : null);
                    }
                }
            }
        }

        return retmap;
    }

    public enum AdditionalInfo
    {
        Gene,
        ANN,
        LOF,
        AA,
        AC,
        AF,
        AN,
        BQ,
        CIGAR,
        DB,
        DP,
        END,
        H2,
        H3,
        MQ,
        MQ0,
        NS,
        SB,
        SOMATIC,
        VALIDATED,
        VE,
        CSQ,
        MAF,
        MAC,
        MA,
        CLIN_risk_factor,
        CLIN_benign,
        CLIN_likely_benign,
        CLIN_uncertain_significance,
        CLIN_histocompatibility,
        CLIN_not_provided,
        CLIN_association,
        CLIN_pathogenic,
        CLIN_likely_pathogenic,
        CLIN_drug_response,
        CLIN_other,
        CLIN_confers_sensitivity,
        CLIN_protective,
        /* dbSNP */
        GENEINFO,
        RS
    }

    public enum Annotation
    {
        Allele, Annotation, Annotation_Impact, Gene_Name, Gene_ID, Feature_Type, Feature_ID, Transcript_BioType, Rank
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
        return this.parseInfo()
                   .get(additionalInfo.name());
    }

    /**
     * Similar to {@link #getInfo()} but returns an {@link Optional}
     * 
     * @param additionalInfo
     * @return
     */
    public Optional<String> getInfoValue(AdditionalInfo additionalInfo)
    {
        return Optional.ofNullable(this.getInfo(additionalInfo));
    }

    /**
     * Returns true, if the given {@link AdditionalInfo} exists in the record. Note: the {@link #getInfo()} function might still return null and the
     * {@link #getInfoValue(AdditionalInfo)} function might still return an empty {@link Optional} for the {@link AdditionalInfo} key, if only the key is
     * present in the record without any value.
     * 
     * @param additionalInfo
     * @return
     */
    public boolean getInfoExists(AdditionalInfo additionalInfo)
    {
        return this.parseInfo()
                   .containsKey(additionalInfo.name());
    }

    /**
     * Similar to {@link #getInfoTokens(AdditionalInfo, char)} with the pipe '|' as separator.
     * 
     * @param additionalInfo
     * @return
     */
    public List<String> getInfoTokens(AdditionalInfo additionalInfo)
    {
        return this.getInfoTokens(additionalInfo, '|');
    }

    /**
     * Returns the n-th value of {@link #getInfoTokens(AdditionalInfo)}. index=0,1,2,...
     * 
     * @param additionalInfo
     * @param valueIndex
     * @return
     */
    public Optional<String> getInfoTokensValue(AdditionalInfo additionalInfo, int valueIndex)
    {
        return ListUtils.optionalFirst(this.getInfoTokens(additionalInfo));
    }

    /**
     * Returns info groups separated by comma (',')
     * 
     * @param additionalInfo
     * @param valueIndex
     * @return
     */
    public Stream<InfoTokenGroup> getInfoTokenGroups(AdditionalInfo additionalInfo)
    {
        return this.getInfoTokenGroups(additionalInfo, ',');
    }

    public Stream<InfoTokenGroup> getInfoTokenGroups(AdditionalInfo additionalInfo, char separator)
    {
        return this.getInfoTokenGroups(additionalInfo, separator, '|');
    }

    public Stream<InfoTokenGroup> getInfoTokenGroups(AdditionalInfo additionalInfo, char groupSeparator, char tokenSeparator)
    {
        return Optional.ofNullable(StringUtils.splitPreserveAllTokens(this.getInfo(additionalInfo), groupSeparator))
                       .map(Arrays::asList)
                       .orElse(Collections.emptyList())
                       .stream()
                       .map(groupInfo ->
                       {
                           List<String> subInfoTokens = this.getInfoTokens(groupInfo, tokenSeparator);
                           return new InfoTokenGroup()
                           {
                               @Override
                               public Optional<String> getValueAt(int index)
                               {
                                   return Optional.ofNullable(subInfoTokens)
                                                  .filter(tokens -> tokens.size() >= index + 1)
                                                  .map(tokens -> tokens.get(index));
                               }
                           };
                       });
    }

    public static interface InfoTokenGroup
    {
        public Optional<String> getValueAt(int index);
    }

    /**
     * Returns the {@link #getInfo()} {@link String} as tokens separated by the given {@link Character}
     * 
     * @param additionalInfo
     * @param separator
     * @return
     */
    public List<String> getInfoTokens(AdditionalInfo additionalInfo, char separator)
    {
        return this.getInfoTokens(this.getInfo(additionalInfo), separator);
    }

    private List<String> getInfoTokens(String info, char separator)
    {
        return Optional.ofNullable(StringUtils.splitPreserveAllTokens(info, separator))
                       .map(Arrays::asList)
                       .orElse(Collections.emptyList());
    }

    /**
     * Returns the gene from the {@link #getInfo(AdditionalInfo)} with the {@link AdditionalInfo#Gene} token if possible, otherwise returns null
     * <br>
     * <br>
     * Example: "BHMT"
     * 
     * @return
     */
    public String getGene()
    {
        String gene = StringUtils.upperCase(org.omnaest.utils.StringUtils.splitToStream(this.getInfo(AdditionalInfo.Gene), "/")
                                                                         .findFirst()
                                                                         .orElse(null));
        if (StringUtils.isBlank(gene))
        {
            gene = StringUtils.upperCase(this.getAnnotation()
                                             .get("Gene_Name"));
        }
        return gene;
    }

    public boolean hasGene()
    {
        return StringUtils.isNotBlank(this.getGene());
    }

    public String getAnnotation(Annotation annotation)
    {
        return this.getAnnotation()
                   .get(annotation.name());
    }

    public Map<String, String> getAnnotation()
    {
        Map<String, String> retmap = Collections.emptyMap();

        String annotation = this.getInfo(AdditionalInfo.ANN);
        if (StringUtils.isNotBlank(annotation))
        {
            List<String> headers = Arrays.asList("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
                                                 "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length",
                                                 "AA.pos / AA.length", "Distance", "ERRORS / WARNINGS / INFO");

            retmap = StreamUtils.merge2(headers.stream(), org.omnaest.utils.StringUtils.splitToStream(annotation, "|"))
                                .filter(element -> element.getFirst() != null)
                                .filter(element -> element.getSecond() != null)
                                .collect(CollectorUtils.toMapByBiElement(() -> new LinkedHashMap<>()));
        }

        return retmap;
    }

    public double getLossOfFunctionPrediction()
    {
        return NumberUtils.toDouble(this.getLossOfFunctionPredictionInfo()
                                        .get("Percent_of_transcripts_affected"));
    }

    public Map<String, String> getLossOfFunctionPredictionInfo()
    {
        Map<String, String> retmap = Collections.emptyMap();

        String lossOfFunctionField = this.getInfo(AdditionalInfo.LOF);
        if (StringUtils.isNotBlank(lossOfFunctionField))
        {
            List<String> headers = Arrays.asList("Gene_Name", "Gene_ID", "Number_of_transcripts_in_gene", "Percent_of_transcripts_affected");

            lossOfFunctionField = lossOfFunctionField.replaceAll("^\\(", "")
                                                     .replaceAll("\\)$", "");
            retmap = StreamUtils.merge2(headers.stream(), org.omnaest.utils.StringUtils.splitToStream(lossOfFunctionField, "|"))
                                .filter(element -> element.getFirst() != null)
                                .filter(element -> element.getSecond() != null)
                                .collect(CollectorUtils.toMapByBiElement(() -> new LinkedHashMap<>()));
        }

        return retmap;
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

    public long getPositionAsLong()
    {
        return NumberUtils.toLong(this.getPosition());
    }

    public boolean hasInsertion()
    {
        int referenceLength = this.getReference()
                                  .length();
        return org.omnaest.utils.StringUtils.splitToStream(this.getAlternativeAlleles(), ",")
                                            .map(String::trim)
                                            .mapToInt(String::length)
                                            .anyMatch(alleleLength -> alleleLength > referenceLength);
    }

    public boolean hasDeletion()
    {
        int referenceLength = this.getReference()
                                  .length();
        return org.omnaest.utils.StringUtils.splitToStream(this.getAlternativeAlleles(), ",")
                                            .map(String::trim)
                                            .mapToInt(String::length)
                                            .anyMatch(alleleLength -> alleleLength < referenceLength);
    }

}
