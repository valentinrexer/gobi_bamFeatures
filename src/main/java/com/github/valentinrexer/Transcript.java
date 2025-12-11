package com.github.valentinrexer;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import augmentedTree.*;
import com.github.valentinrexer.utils.BamFeatureUtils;

public class Transcript {
    private final String transcriptId;
    private final String geneId;
    private final char strand;
    private final List<Exon> exons;
    private List<Region> exonVector;
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

    private void computeExonVector() {
        var exonVector = new ArrayList<Region>();
        for (Exon exon : exons)
            exonVector.add(exon.toRegion());

        BamFeatureUtils.mergeVector(exonVector);
        this.exonVector = exonVector;
    }

    public List<Region> getExonVector() {
        if (exonVector == null) computeExonVector();
        return exonVector;
    }

    public List<Region> getExonRegionsForInterval(Region interval) {
        List<Region> regions = new ArrayList<>();

        for (Exon exon : exons) {
            if(exon.toRegion().intersects(interval))
                regions.add(new Region(Math.max(interval.start(), exon.getStart()),
                        Math.min(interval.end(), exon.getEnd()))
                );
        }

        return regions;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Transcript that = (Transcript) o;

        if (strand != that.strand) return false;
        if (!transcriptId.equals(that.transcriptId)) return false;
        if (!geneId.equals(that.geneId)) return false;

        // Compare exons structurally: same count + same coordinates
        if (this.exons.size() != that.exons.size()) return false;

        // Ensure consistent ordering before structural comparison
        List<Exon> ex1 = new ArrayList<>(this.exons);
        List<Exon> ex2 = new ArrayList<>(that.exons);

        ex1.sort(Comparator.comparingInt(Exon::getStart));
        ex2.sort(Comparator.comparingInt(Exon::getStart));

        for (int i = 0; i < ex1.size(); i++) {
            Exon a = ex1.get(i);
            Exon b = ex2.get(i);

            if (a.getStart() != b.getStart() || a.getEnd() != b.getEnd())
                return false;
        }

        return true;
    }

    @Override
    public int hashCode() {
        int result = transcriptId.hashCode();
        result = 31 * result + geneId.hashCode();
        result = 31 * result + Character.hashCode(strand);

        // Include exon structure in hash
        List<Exon> sorted = new ArrayList<>(exons);
        sorted.sort(Comparator.comparingInt(Exon::getStart));

        for (Exon e : sorted) {
            result = 31 * result + e.getStart();
            result = 31 * result + e.getEnd();
        }

        return result;
    }

    @Override
    public String toString() {
        return "Transcript{id=" + transcriptId +
                ", gene=" + geneId +
                ", strand=" + strand +
                ", exons=" + exons.size() + "}";
    }
}

