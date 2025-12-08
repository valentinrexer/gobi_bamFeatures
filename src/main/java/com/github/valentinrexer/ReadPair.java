package com.github.valentinrexer;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.*;
import java.util.stream.Collectors;

public class ReadPair {
    private final SAMRecord firstRecord;
    private final SAMRecord lastRecord;
    private final Region firstRecordRegion;
    private final Region lastRecordRegion;
    private final List<Region> regionVector;

    public ReadPair(SAMRecord firstRecord, SAMRecord lastRecord) {
        this.firstRecord = firstRecord;
        this.lastRecord = lastRecord;

        List<Region> regionVector = getRegionVector(firstRecord, false);
        regionVector.addAll(getRegionVector(lastRecord, false));
        regionVector = mergeVector(regionVector);
        this.regionVector = regionVector;

        firstRecordRegion = new Region(firstRecord.getAlignmentStart(), firstRecord.getAlignmentEnd());
        lastRecordRegion = new Region(lastRecord.getAlignmentStart(), lastRecord.getAlignmentEnd());
    }

    public String process(GtfData gtfData, Boolean frStrand) {
        Integer nSplit = getNSplit();
        if (nSplit == null) return firstRecord.getReadName() + "\tsplit-inconsistent:true";

        int mm = getMismatches();
        int clipped = getTotalClipped();
        var geneLvl = getGeneAnnotation(gtfData, frStrand);

        int gCount;
        String outString;
        if (geneLvl.getFirst().level() == GenicLevel.INTERGENIC) {
            gCount = 0;
            var geneDistance = "";
            outString = "";
        }
        else {
            gCount = geneLvl.size();
            StringBuilder associatedGenesString = new StringBuilder();
            for (GenicLevelContainer container : geneLvl)
                associatedGenesString.append(container.annotationString()).append("|");

            outString = associatedGenesString.substring(0, associatedGenesString.length() - 1);
        }
        System.out.println("gcount:" + gCount + "\t" + outString);

        return "";
    }
    
    private int getGeneDistance(String chr, GtfData gtfData, Boolean frStrand) {
        var readBounds = getMinStartMaxEnd();
        var leftNeighbor = gtfData.getLeftNeighbor(chr, readBounds.start(), readBounds.end(), frStrand);
        var rightNeighbor = gtfData.getRightNeighbor(chr, readBounds.start(), readBounds.end(), frStrand);

        int minLeftDist = Integer.MAX_VALUE, minRightDist = Integer.MAX_VALUE;
        for (Gene neighbor : leftNeighbor) {
            var dist = readBounds.getStart() - neighbor.getEnd();
            if (dist < 0) return 0;

            minRightDist = Math.min(minRightDist, dist);
        }

        for (Gene neighbor : rightNeighbor) {
            var dist = readBounds.end() - neighbor.getStart();
            if (dist < 0) return 0;
            minLeftDist = Math.min(minLeftDist, dist);
        }

        return Math.min(minLeftDist, minRightDist);
    } 

    private List<GenicLevelContainer> getGeneAnnotation(GtfData gtfData, Boolean frStrand) {
        var candidates = getCandidateGenes(gtfData, frStrand);
        var genicLevelMapping = new HashMap<GenicLevel, List<GenicLevelContainer>>();

        for (Gene candidate : candidates) {
            GenicLevelContainer container = getGenicLevel(candidate);

            genicLevelMapping
                    .computeIfAbsent(container.level(), k -> new ArrayList<>())
                    .add(container);
        }

        if (genicLevelMapping.containsKey(GenicLevel.TRANSCRIPTOMIC)) return genicLevelMapping.get(GenicLevel.TRANSCRIPTOMIC);
        if (genicLevelMapping.containsKey(GenicLevel.MERGED_TRANSCRIPTOMIC)) return genicLevelMapping.get(GenicLevel.MERGED_TRANSCRIPTOMIC);
        if (genicLevelMapping.containsKey(GenicLevel.INTRONIC)) return genicLevelMapping.get(GenicLevel.INTRONIC);
        return List.of(new GenicLevelContainer(GenicLevel.INTERGENIC, null));
    }

    private GenicLevelContainer getGenicLevel(Gene candidateGene) {
        List<Transcript> matchingTranscripts = new ArrayList<>();

        for (Transcript transcript : candidateGene.getTranscripts()) {
            var tree = transcript.getExonIntervalTree();
            var covers = true;

            for (Region vector : regionVector) {
                var exonsCoveringVector = tree.getIntervalsSpanning(vector.start(), vector.end(), new ArrayList<>());
                if (exonsCoveringVector.isEmpty()) {
                    covers = false;
                    break;
                }
            }
            if (!covers) continue;
            matchingTranscripts.add(transcript);
        }

        if (!matchingTranscripts.isEmpty()) {
            StringBuilder ret = new StringBuilder(candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":");
            for (Transcript transcript : matchingTranscripts)
                ret.append(transcript.getGeneId()).append(",");

            return new GenicLevelContainer(GenicLevel.TRANSCRIPTOMIC, ret.substring(0, ret.length() - 1));
        }

        var tree = candidateGene.getMergedTranscriptTree();
        for (Region vector : regionVector) {
            var exonsCoveringVector = tree.getIntervalsSpanning(vector.start(), vector.end(), new ArrayList<>());
            if (exonsCoveringVector.isEmpty())
                return new GenicLevelContainer(GenicLevel.INTRONIC, candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":INTRON");
        }

        return new GenicLevelContainer(GenicLevel.MERGED_TRANSCRIPTOMIC, candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":MERGED");
    }

    private ReadPairType determineReadPairTypeForSingleGene(Gene gene) {
        Region interval = getMinStartMaxEnd();
        if (interval.start() >= gene.getStart() && interval.end() < gene.getEnd()) return ReadPairType.GENIC;
        else if (interval.start() < gene.getStart() && interval.end() > gene.getEnd()) return ReadPairType.COVERS_WHOLE_GENE;
        else return ReadPairType.INTERGENIC;
    }

    private List<Gene> getCandidateGenes(GtfData gtfData, Boolean frStrand) {
        String chr = firstRecord.getReferenceName();

        var readBounds = getMinStartMaxEnd();

        HashSet<Gene> candidates = new HashSet<>(gtfData.getContainingGenes(chr, readBounds.start(), readBounds.end(), frStrand));
        
        return new ArrayList<>(candidates);
    }

    private static List<Region> mergeVector(List<Region> vector) {
        if (vector.isEmpty()) return List.of();

        vector.sort(Comparator.comparingInt(Region::start)
                .thenComparingInt(Region::end));

        var merged = new ArrayList<Region>();
        Region prev = vector.get(0);

        for (int i = 1; i < vector.size(); i++) {
            Region next = vector.get(i);

            if (prev.end() >= next.start() - 1) {
                prev = new Region(prev.start(), Math.max(prev.end(), next.end()));
            } else {
                merged.add(prev);
                prev = next;
            }
        }
        merged.add(prev);

        return merged;
    }

    private static List<Region> getRegionVector(SAMRecord record, Boolean intron) {
        List<Region> regions = new ArrayList<>();

        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        var pointer = record.getAlignmentStart();

        for (CigarElement cigarElement : cigarElements) {
            var length = cigarElement.getLength();
            CigarOperator operator = cigarElement.getOperator();

            switch (operator) {
                case M:
                case EQ:
                case X:
                case D:
                    if (intron == null || !intron) regions.add(new Region(pointer, pointer + length - 1));
                    pointer += length;
                    break;

                case N:
                    if (intron == null || intron) regions.add(new Region(pointer, pointer + length - 1));
                    pointer += length;
                    break;

                case I:
                case S:
                case H:
                case P:
                    break;

                default:
                    throw new IllegalStateException("Unknown operator: " + operator);
            }
        }

        return regions;
    }

    private Integer getNSplit() {
        Region firstRecordRegion = new Region(firstRecord.getAlignmentStart(), firstRecord.getAlignmentEnd());
        Region lastRecordRegion = new Region(lastRecord.getAlignmentStart(), lastRecord.getAlignmentEnd());
        Region intersectingRegion = firstRecordRegion.getIntersectingRegion(lastRecordRegion);

        if (!(intersectingRegion == null))
            if (!isConsistentSplit()) return null;


        int nCount = 0;
        for (SAMRecord record : List.of(firstRecord, lastRecord)){
            for (CigarElement element : record.getCigar().getCigarElements())
                if (element.getOperator() == CigarOperator.N) nCount++;
        }
        return nCount;
    }

    private boolean isConsistentSplit() {
        Region overlap = firstRecordRegion.getIntersectingRegion(lastRecordRegion);

        // No shared region → they cannot contradict each other → treat as consistent
        if (overlap == null) return true;

        // Get introns as Sets<Region>
        Set<Region> firstIntrons = getIntrons(firstRecord);
        Set<Region> lastIntrons  = getIntrons(lastRecord);

        // Keep only introns inside the shared region
        firstIntrons = firstIntrons.stream()
                .filter(i -> i.intersects(overlap))
                .collect(Collectors.toSet());

        lastIntrons = lastIntrons.stream()
                .filter(i -> i.intersects(overlap))
                .collect(Collectors.toSet());

        // Consistent = both sets identical
        return firstIntrons.equals(lastIntrons);
    }


    private static HashSet<Region> getIntrons(SAMRecord record) {
        var introns = new HashSet<Region>();
        List<AlignmentBlock> blocks = record.getAlignmentBlocks();

        if (blocks.size() < 2) return introns;

        for (int i = 0; i < blocks.size() - 1; i++) {
            AlignmentBlock a = blocks.get(i);
            AlignmentBlock b = blocks.get(i + 1);

            int intronStart = a.getReferenceStart() + a.getLength();
            int intronEnd   = b.getReferenceStart() - 1;

            if (intronEnd >= intronStart) {
                introns.add(new Region(intronStart, intronEnd));
            }
        }

        return introns;
    }


    private int getMismatches() {
        Integer nm1 = firstRecord.getIntegerAttribute("NM");
        Integer nm2 = lastRecord.getIntegerAttribute("NM");

        if (nm1 == null) nm1 = 0;
        if (nm2 == null) nm2 = 0;

        return nm1 + nm2;
    }

    private int getTotalClipped() {
        return getClippingCount(firstRecord) + getClippingCount(lastRecord);
    }

    private static int getClippingCount(SAMRecord record) {
        List<CigarElement> cigar = record.getCigar().getCigarElements();
        int clipped = 0;

        for (CigarElement element : cigar) {
            CigarOperator operator = element.getOperator();

            if (operator == CigarOperator.S || operator == CigarOperator.H) {
                clipped += element.getLength();
            }
        }

        return clipped;
    }

    private Region getMinStartMaxEnd() {
        int minPos = Integer.MAX_VALUE;
        int maxPos = Integer.MIN_VALUE;

        for (Region region : regionVector) {
            minPos = Math.min(minPos, region.start());
            maxPos = Math.max(maxPos, region.end());
        }

        return new Region(minPos, maxPos);
    }

    @Override
    public String toString() {
        return firstRecord.getReadName() + ", " + lastRecord.getReadName();
    }
}
