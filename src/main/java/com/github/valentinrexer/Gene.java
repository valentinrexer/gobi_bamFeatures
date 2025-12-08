package com.github.valentinrexer;

import augmentedTree.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Gene implements Interval {
    private final String geneId;
    private final String geneName;
    private final String geneBiotype;
    private final char strand;
    private final String chromosome;
    private final List<Transcript> transcripts;
    private IntervalTree<Region> mergedTranscriptTree;
    private int start = Integer.MAX_VALUE;
    private int end = Integer.MIN_VALUE;

    public Gene(String geneId, String geneName, String geneBiotype, char strand, String chromosome) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.geneBiotype = geneBiotype;
        this.strand = strand;
        this.chromosome = chromosome;
        this.transcripts = new ArrayList<>();
    }

    public void addTranscript(Transcript transcript) {
        if(!transcripts.contains(transcript))
            transcripts.add(transcript);
    }

    public void computeBoundaries() {
        for (Transcript transcript : transcripts) {
            transcript.computeBoundaries();
            start = Math.min(start, transcript.getStart());
            end = Math.max(end, transcript.getEnd());
        }
    }

    public void sortExonsForEachTranscript() {
        transcripts.forEach(Transcript::sortExons);
    }

    public String getGeneId() { return geneId; }
    public String getGeneName() { return geneName; }
    public String getGeneBiotype() { return geneBiotype; }
    public char getStrand() { return strand; }
    public String getChromosome() { return chromosome; }

    public List<Transcript> getTranscripts() {
        return Collections.unmodifiableList(transcripts);
    }

    public IntervalTree<Region> getMergedTranscriptTree() {
        if (mergedTranscriptTree == null)
            computeMergedTranscriptTree();

        return mergedTranscriptTree;
    }

    private void computeMergedTranscriptTree() {
        mergedTranscriptTree = new IntervalTree<>();

        for  (Transcript transcript : transcripts) {
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

