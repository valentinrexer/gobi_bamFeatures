package com.github.valentinrexer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import augmentedTree.*;

public class Transcript {
    private final String transcriptId;
    private final String geneId;
    private final char strand;
    private final List<Exon> exons;
    private int start = Integer.MAX_VALUE;
    private int end = Integer.MIN_VALUE;

    public Transcript(String transcriptId, String geneId, char strand) {
        this.transcriptId = transcriptId;
        this.geneId = geneId;
        this.strand = strand;
        this.exons = new ArrayList<>();
    }

    public void addExon(Exon exon) {
        exons.add(exon);
    }

    public void computeBoundaries() {
        for (Exon exon : exons) {
            start = Math.min(start, exon.getStart());
            end = Math.max(end, exon.getEnd());
        }
    }

    public void sortExons() {
        exons.sort((e1, e2) -> Integer.compare(e1.getStart(), e2.getStart()));
    }

    public String getTranscriptId() { return transcriptId; }
    public String getGeneId() { return geneId; }
    public char getStrand() { return strand; }

    public List<Exon> getExons() {
        return Collections.unmodifiableList(exons);
    }

    public int getStart() { return start; }
    public int getEnd() { return end; }

    public IntervalTree<Region> getExonIntervalTree() {
        var tree = new IntervalTree<Region>();

        for (Exon exon : exons)
            tree.add(new  Region(exon.getStart(), exon.getEnd()));

        return tree;
    }

    @Override
    public String toString() {
        return "Transcript{id=" + transcriptId +
                ", gene=" + geneId +
                ", strand=" + strand +
                ", exons=" + exons.size() + "}";
    }
}

