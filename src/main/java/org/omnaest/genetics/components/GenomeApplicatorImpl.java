package org.omnaest.genetics.components;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.omnaest.genetics.domain.VCFData.GenomeApplicator;
import org.omnaest.genetics.domain.VCFData.Replacements;
import org.omnaest.genetics.domain.VCFRecord;
import org.omnaest.genetics.domain.VCFRecord.SampleFields;
import org.omnaest.genetics.domain.VCFRecord.SampleFields.GenoType;
import org.omnaest.genetics.translator.domain.CodeAndPosition;
import org.omnaest.genetics.translator.domain.NucleicAcidCode;
import org.omnaest.utils.ConsumerUtils;
import org.omnaest.utils.SetUtils;
import org.omnaest.utils.element.lar.UnaryLeftAndRight;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GenomeApplicatorImpl implements GenomeApplicator
{
    private static final Logger LOG = LoggerFactory.getLogger(GenomeApplicatorImpl.class);

    private final Map<String, List<VCFRecord>> chromosomeToRecords;

    public static class ReplacementsImpl implements Replacements
    {
        private Map<Integer, Set<UnaryLeftAndRight<NucleicAcidCode>>> alleleToReplacements = new HashMap<>();

        public void addReplacementForAllele(int allele, UnaryLeftAndRight<NucleicAcidCode> replacement)
        {
            this.alleleToReplacements.computeIfAbsent(allele, a -> new HashSet<>())
                                     .add(replacement);
        }

        @Override
        public Set<UnaryLeftAndRight<NucleicAcidCode>> getReplacementForAllele(int allele)
        {
            return this.alleleToReplacements.get(allele);
        }

        @Override
        public int getMaxAlleleIndex()
        {
            return this.alleleToReplacements.keySet()
                                            .stream()
                                            .mapToInt(v -> v)
                                            .max()
                                            .orElse(-1);
        }

        @Override
        public boolean hasReplacementForAllele(int allele)
        {
            return this.alleleToReplacements.containsKey(allele) && !this.alleleToReplacements.get(allele)
                                                                                              .isEmpty();
        }

    }

    public GenomeApplicatorImpl(Map<String, List<VCFRecord>> chromosomeToRecords)
    {
        this.chromosomeToRecords = chromosomeToRecords;
    }

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
    public Map<Long, Replacements> getPositionToReplacementForChromosome(String chromosome)
    {
        return this.determinePositionToReplacement(this.chromosomeToRecords.getOrDefault(StringUtils.upperCase(chromosome), Collections.emptyList()));
    }

    @Override
    public Stream<ChromosomeAndPositionReplacement> getPositionToReplacements()
    {
        return this.chromosomeToRecords.keySet()
                                       .stream()
                                       .map(chromosome -> new ChromosomeAndPositionReplacement()
                                       {
                                           @Override
                                           public Map<Long, Replacements> getPositionToReplacement()
                                           {
                                               return GenomeApplicatorImpl.this.getPositionToReplacementForChromosome(chromosome);
                                           }

                                           @Override
                                           public String getChromosome()
                                           {
                                               return chromosome;
                                           }
                                       });
    }

    private Map<Long, Replacements> determinePositionToReplacement(List<VCFRecord> records)
    {
        Map<Long, Replacements> positionToReplacements = new ConcurrentHashMap<>();

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
                    NucleicAcidCode right = ii < alternativeAlleles.length() ? NucleicAcidCode.valueOf(alternativeAlleles.charAt(ii)) : null;
                    ReplacementsImpl replacements = (ReplacementsImpl) positionToReplacements.computeIfAbsent(currentPosition, c -> new ReplacementsImpl());

                    SampleFields sampleFields = vcfRecord.parseSampleFields();
                    if (sampleFields.hasGenoType(GenoType.ALTERNATIVE_BOTH))
                    {
                        replacements.addReplacementForAllele(0, new UnaryLeftAndRight<NucleicAcidCode>(left, right));
                        replacements.addReplacementForAllele(1, new UnaryLeftAndRight<NucleicAcidCode>(left, right));
                    }
                    else if (sampleFields.hasGenoType(GenoType.REFERENCE_AND_ALTERNATIVE))
                    {
                        if (sampleFields.hasGenoType(GenoType.SECOND_ALLELE_ALTERNATIVE))
                        {
                            replacements.addReplacementForAllele(1, new UnaryLeftAndRight<NucleicAcidCode>(left, right));
                        }
                        else
                        {
                            replacements.addReplacementForAllele(0, new UnaryLeftAndRight<NucleicAcidCode>(left, right));
                        }
                    }
                    else
                    {
                        int allele = Math.min(1, replacements.getMaxAlleleIndex() + 1);
                        replacements.addReplacementForAllele(allele, new UnaryLeftAndRight<NucleicAcidCode>(left, right));
                    }
                }
            }
        }

        return positionToReplacements;
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
                return this.applyToChromosomeCodeAndPositionSequence(chromosome, sequence.map(code -> new CodeAndPosition<>(code, position.getAndIncrement())))
                           .map(cap -> cap.getCode());
            }

            @Override
            public Stream<CodeAndPosition<NucleicAcidCode>> applyToChromosomeCodeAndPositionSequence(String chromosome,
                                                                                                     Stream<CodeAndPosition<NucleicAcidCode>> sequence)
            {
                Map<Long, Replacements> positionToReplacement = GenomeApplicatorImpl.this.getPositionToReplacementForChromosome(chromosome);

                AtomicLong position = new AtomicLong(-1);
                return sequence.peek(ConsumerUtils.consumeOnce(code ->
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
                                   Replacements replacementsHolder = positionToReplacement.get(currentPosition);
                                   Set<UnaryLeftAndRight<NucleicAcidCode>> replacements = replacementsHolder != null
                                           ? replacementsHolder.getReplacementForAllele(allele)
                                           : null;
                                   if (replacements != null)
                                   {
                                       UnaryLeftAndRight<NucleicAcidCode> replacement = SetUtils.first(replacements);

                                       if (replacements.size() > 1)
                                       {
                                           LOG.warn("More than one replacement for chromosome position and allele available: " + chromosome + ":" + position
                                                   + " ( allele " + allele + " )");
                                           LOG.warn(replacements.toString());
                                       }

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
                    throw new IllegalStateException("Reference code did not match: " + code + "<->" + referenceCode + " at position: " + currentPosition);
                }
            }
        };
    }

    @Override
    public int getNumberOfAlleles()
    {
        return 1 + this.getPositionToReplacements()
                       .flatMap(capr -> capr.getPositionToReplacement()
                                            .values()
                                            .stream()
                                            .map(replacements -> replacements.getMaxAlleleIndex()))
                       .mapToInt(v -> v)
                       .max()
                       .getAsInt();
    }
}