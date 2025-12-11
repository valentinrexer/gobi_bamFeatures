package com.github.valentinrexer;

import augmentedTree.*;
import com.github.valentinrexer.utils.BamFeatureUtils;

import java.util.*;

public class Gene implements Interval {
    private final String geneId;
    private final String geneName;
    private final String geneBiotype;
    private final char strand;
    private final String chromosome;
    private final HashMap<String, Transcript> transcripts;
    private IntervalTree<Region> mergedTranscriptome;
    private int start = Integer.MAX_VALUE;
    private int end = Integer.MIN_VALUE;

    public Gene(String geneId, String geneName, String geneBiotype, char strand, String chromosome) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.geneBiotype = geneBiotype;
        this.strand = strand;
        this.chromosome = chromosome;
        this.transcripts = new HashMap<>();
    }

    public void addTranscript(Transcript transcript) {
        transcripts.put(transcript.getTranscriptId(), transcript);
    }

    public Transcript getTranscript(String transcriptId) {
        if (!transcripts.containsKey(transcriptId)) return null;
        return transcripts.get(transcriptId);
    }

    public void computeBoundaries() {
        for (Transcript transcript : transcripts.values()) {
            transcript.computeBoundaries();
            start = Math.min(start, transcript.getStart());
            end = Math.max(end, transcript.getEnd());
        }
    }

    public void sortExonsForEachTranscript() {
        transcripts.values().forEach(Transcript::sortExons);
    }

    public String getGeneId() { return geneId; }
    public String getGeneName() { return geneName; }
    public String getGeneBiotype() { return geneBiotype; }
    public char getStrand() { return strand; }
    public String getChromosome() { return chromosome; }

    public List<Transcript> getTranscripts() {
        return transcripts.values().stream().toList();
    }

    public void computeMergedTranscriptome() {
        List<Region> regions = new ArrayList<>();

        for (Transcript transcript : transcripts.values()) {
            regions.addAll(transcript.getExonVector());
        }

        regions = BamFeatureUtils.mergeVector(regions);

        mergedTranscriptome = new IntervalTree<>();
        mergedTranscriptome.addAll(regions);
    }

    public List<Region> getMergedTranscriptomeForInterval(Region interval) {
        if (mergedTranscriptome == null) computeMergedTranscriptome();

        List<Region> intersectingRegions = mergedTranscriptome.getIntervalsIntersecting(interval.start(), interval.end(), new ArrayList<>());
        List<Region> trimmedRegions = new ArrayList<>();

        for (Region intersectingRegion : intersectingRegions) {
            trimmedRegions.add(new Region(Math.max(interval.start(), intersectingRegion.start()),
                    Math.min(interval.end(), intersectingRegion.end())));
        }
        return trimmedRegions;
    }

    @Override
    public int getStart() { return start; }

    @Override
    public int getStop() { return end; }

    public int getEnd() { return end; }

    @Override
    public String toString() {
        return "Gene{id=" + geneId +
                ", name=" + geneName +
                ", biotype=" + geneBiotype +
                ", strand=" + strand +
                ", transcripts=" + transcripts.size() +
                ", start=" +  start +
                ", end=" + end +  "}";
    }
}

