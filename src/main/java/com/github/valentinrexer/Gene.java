package com.github.valentinrexer;

import augmentedTree.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Gene implements Interval {
    private final String geneId;
    private final String geneName;
    private final String geneBiotype;
    private final char strand;
    private final String chromosome;
    private final HashMap<String, Transcript> transcripts;
    private IntervalTree<Region> mergedTranscriptTree;
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
        return Collections.unmodifiableList(transcripts.values().stream().toList());
    }

    public IntervalTree<Region> getMergedTranscriptTree() {
        if (mergedTranscriptTree == null)
            computeMergedTranscriptTree();

        return mergedTranscriptTree;
    }

    private void computeMergedTranscriptTree() {
        mergedTranscriptTree = new IntervalTree<>();

        for  (Transcript transcript : transcripts.values()) {
            for (Exon exon : transcript.getExons()) {
                mergedTranscriptTree.add(new Region(exon.getStart(), exon.getEnd()));
            }
        }
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

