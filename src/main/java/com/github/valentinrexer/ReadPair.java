package com.github.valentinrexer;

import com.github.valentinrexer.utils.BamFeatureUtils;
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
    private final List<Region> pairRegionVector;
    private final List<Region> regionVectorFirst;
    private final List<Region> regionVectorLast;
    private final Set<Region> intronsFirst;
    private final Set<Region> intronsLast;
    private final String chromosome;
    private final boolean strand;

    public ReadPair(SAMRecord firstRecord, SAMRecord lastRecord) {
        this.firstRecord = firstRecord;
        this.lastRecord = lastRecord;

        regionVectorFirst = getRegionVector(firstRecord);
        regionVectorLast = getRegionVector(lastRecord);
        strand = !firstRecord.getReadNegativeStrandFlag();

        List<Region> regionVector = new ArrayList<>(regionVectorFirst);
        regionVector.addAll(regionVectorLast);
        regionVector = BamFeatureUtils.mergeVector(regionVector);
        this.pairRegionVector = regionVector;
        this.chromosome = firstRecord.getReferenceName();

        firstRecordRegion = new Region(firstRecord.getAlignmentStart(), firstRecord.getAlignmentEnd());
        lastRecordRegion = new Region(lastRecord.getAlignmentStart(), lastRecord.getAlignmentEnd());

        this.intronsFirst = getIntrons(firstRecord);
        this.intronsLast = getIntrons(lastRecord);
    }

    public String process(TreeGtf treeGtf, Boolean frStrand, PcrIndexMap pcrIndexMap) {
        Integer nSplit = getNSplit();
        if (nSplit == null) return firstRecord.getReadName() + "\tsplit-inconsistent:true";

        int mm = getMismatches();
        int clipping = getTotalClipped();
        var geneLvl = getGeneAnnotation(treeGtf, frStrand);
        int gCount;
        String geneOutputString;

        if (geneLvl.getFirst().level() == GenicLevel.INTERGENIC) {
            gCount = 0;
            var geneDistance = getGeneDistance(chromosome, treeGtf, frStrand);
            geneOutputString = "gdist:" + geneDistance;

            var hasAntisenseGene = hasAntiSenseGene(treeGtf, frStrand);
            geneOutputString += "\tantisense:" + hasAntisenseGene;
        }
        else {
            gCount = geneLvl.size();
            StringBuilder associatedGenesString = new StringBuilder();
            for (GenicLevelContainer container : geneLvl)
                associatedGenesString.append(container.annotationString()).append("|");

            geneOutputString = associatedGenesString.substring(0, associatedGenesString.length() - 1);
        }

        Boolean indexStrand = frStrand == null ? null : (frStrand == strand);
        int pcrIndex = pcrIndexMap.getPcrIndex(pairRegionVector, indexStrand);

        return firstRecord.getReadName() +
                "\tmm:" + mm +
                "\tclipping:" + clipping +
                "\tgcount:"+ gCount +
                "\tnsplit:" + nSplit +
                "\t" + geneOutputString +
                "\tpcrindex:" + pcrIndex;
    }

    private boolean hasAntiSenseGene(TreeGtf treeGtf, Boolean frStrand) {
        if (frStrand == null) return false;
        var annotation = getGeneAnnotation(treeGtf, !frStrand);
        return annotation.getFirst().level() != GenicLevel.INTERGENIC;
    }
    
    private int getGeneDistance(String chr, TreeGtf treeGtf, Boolean frStrand) {
        var readBounds = getMinStartMaxEnd();
        var leftNeighbor = treeGtf.getLeftNeighbor(chr, readBounds.start(), readBounds.end(), frStrand);
        var rightNeighbor = treeGtf.getRightNeighbor(chr, readBounds.start(), readBounds.end(), frStrand);

        int minLeftDist = Integer.MAX_VALUE, minRightDist = Integer.MAX_VALUE;
        for (Gene neighbor : leftNeighbor) {
            var dist = readBounds.getStart() - neighbor.getEnd();
            if (dist < 0) return 0;

            minRightDist = Math.min(minRightDist, dist);
        }

        for (Gene neighbor : rightNeighbor) {
            var dist = neighbor.getStart() - readBounds.end();
            if (dist < 0) return 0;
            minLeftDist = Math.min(minLeftDist, dist);
        }

        return Math.min(minLeftDist, minRightDist);
    } 

    private List<GenicLevelContainer> getGeneAnnotation(TreeGtf treeGtf, Boolean frStrand) {
        Boolean lookupStrand =
                frStrand == null ? null : (frStrand == strand);

        var candidates = getCandidateGenes(treeGtf, lookupStrand);
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
        String transcriptomicString = isTranscriptomic(candidateGene);
        if (transcriptomicString != null) return new GenicLevelContainer(GenicLevel.TRANSCRIPTOMIC, transcriptomicString);

        String mergedTranscriptomicString = isMergedTranscriptomic(candidateGene);
        if (mergedTranscriptomicString != null) return new GenicLevelContainer(GenicLevel.MERGED_TRANSCRIPTOMIC, mergedTranscriptomicString);

        else return new GenicLevelContainer(GenicLevel.INTRONIC, candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":INTRON");
    }

    private String isTranscriptomic(Gene candidateGene) {
        List<Transcript> matchingTranscripts = new ArrayList<>();

        var exonVectorFirstRead = new HashSet<>(regionVectorFirst);
        var exonVectorLastRead = new HashSet<>(regionVectorLast);

        for (Transcript transcript : candidateGene.getTranscripts()) {
            var exonVectorFirstInterval = new HashSet<>(transcript.getExonRegionsForInterval(firstRecordRegion));
            var exonVectorLastInterval = new HashSet<>(transcript.getExonRegionsForInterval(lastRecordRegion));

            if (exonVectorFirstInterval.equals(exonVectorFirstRead) && exonVectorLastInterval.equals(exonVectorLastRead)) {
                matchingTranscripts.add(transcript);
            }
        }

        if (matchingTranscripts.isEmpty()) return null;

        StringBuilder transcriptomicString = new StringBuilder(candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":");
        for (Transcript transcript : matchingTranscripts)
            transcriptomicString.append(transcript.getTranscriptId()).append(",");

        return transcriptomicString.substring(0, transcriptomicString.length() - 1);
    }

    private String isMergedTranscriptomic(Gene candidateGene) {
        for (Region block : regionVectorFirst) {
            var mergedTranscriptBlock = candidateGene.getMergedTranscriptomeForInterval(block);
            if (mergedTranscriptBlock.size() != 1) return null;
            if (!mergedTranscriptBlock.getFirst().equals(block)) return null;
        }

        for (Region block : regionVectorLast) {
            var mergedTranscriptBlock = candidateGene.getMergedTranscriptomeForInterval(block);
            if (mergedTranscriptBlock.size() != 1) return null;
            if (!mergedTranscriptBlock.getFirst().equals(block)) return null;
        }

        return candidateGene.getGeneId() + "," + candidateGene.getGeneBiotype() + ":MERGED";
    }

    private List<Gene> getCandidateGenes(TreeGtf treeGtf, Boolean strand) {
        String chr = firstRecord.getReferenceName();

        var readBounds = getMinStartMaxEnd();
        HashSet<Gene> candidates = new HashSet<>(treeGtf.getContainingGenes(chr, readBounds.start(), readBounds.end(), strand));

        return new ArrayList<>(candidates);
    }

    private Integer getNSplit() {
        Region overlap = firstRecordRegion.getIntersectingRegion(lastRecordRegion);

        Set<Region> firstIntrons = new HashSet<>(intronsFirst);
        Set<Region> lastIntrons  = new HashSet<>(intronsLast);

        if (overlap != null) {
            var firstIntronsOverlap = firstIntrons.stream()
                    .filter(i -> i.intersects(overlap))
                    .collect(Collectors.toSet());

            var lastIntronsOverlap = lastIntrons.stream()
                    .filter(i -> i.intersects(overlap))
                    .collect(Collectors.toSet());

            if (!firstIntronsOverlap.equals(lastIntronsOverlap))
                return null;
        }

        firstIntrons.addAll(lastIntrons);
        return firstIntrons.size();
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

        for (Region region : pairRegionVector) {
            minPos = Math.min(minPos, region.start());
            maxPos = Math.max(maxPos, region.end());
        }

        return new Region(minPos, maxPos);
    }

    private static HashSet<Region> getIntrons(SAMRecord record) {
        var introns = new HashSet<Region>();
        List<AlignmentBlock> blocks = record.getAlignmentBlocks();

        if (blocks.size() < 2) return introns;

        for (int i = 0; i < blocks.size() - 1; i++) {
            AlignmentBlock a = blocks.get(i);
            AlignmentBlock b = blocks.get(i + 1);

            int intronStart = a.getReferenceStart() + a.getLength();
            int intronEnd = b.getReferenceStart() - 1;

            if (intronEnd >= intronStart) {
                introns.add(new Region(intronStart, intronEnd));
            }
        }

        return introns;
    }

    private static List<Region> getRegionVector(SAMRecord record) {
        List<Region> regions = new ArrayList<>();

        for (AlignmentBlock block : record.getAlignmentBlocks()) {
            regions.add(new Region(block.getReferenceStart(), block.getReferenceStart() + block.getLength() - 1));
        }

        return BamFeatureUtils.mergeVector(regions);
    }

    @Override
    public String toString() {
        return firstRecord.getReadName() + ", " + lastRecord.getReadName();
    }
}
