package com.github.valentinrexer;

import augmentedTree.IntervalTree;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class TreeGtf {
    private final HashMap<String, Gene> genes = new HashMap<>();
    private final Map<Boolean, Map<String, IntervalTree<Gene>>> geneTrees = new HashMap<>();

    public void readInGffFile(Path filePath, Boolean frStrand) {
        try (BufferedReader br = Files.newBufferedReader(filePath)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#') continue;

                GffLine gffLine = new GffLine(line);
                String gene_id = gffLine.getAttribute("gene_id");
                if (gene_id == null) continue;

                String transcript_id = gffLine.getAttribute("transcript_id");
                char strand = gffLine.getStrand();

                if (!genes.containsKey(gene_id)) {
                    genes.put(gene_id, new Gene(
                            gene_id,
                            gffLine.getAttribute("gene_name"),
                            gffLine.getAttribute("gene_biotype"),
                            strand,
                            gffLine.getSeqId()
                    ));
                } else {
                    if (genes.get(gene_id).getTranscript(transcript_id) == null) {
                        genes.get(gene_id).addTranscript(new Transcript(
                                transcript_id,
                                gene_id,
                                strand
                        ));
                    }

                    if (gffLine.getType().equals("exon")) {
                        genes.get(gene_id).getTranscript(transcript_id).addExon(
                                new Exon(gffLine.getStart(), gffLine.getEnd(), transcript_id)
                        );
                    }
                }
            }

            genes.values().forEach(Gene::computeBoundaries);
            buildIntervalTrees(frStrand);
        } catch (RuntimeException | IOException e) {
            System.err.println("Error occurred during initialization: " + e.getMessage());
        }
    }

    public void buildIntervalTrees(Boolean frStrand) {
        genes.values().forEach(Gene::computeBoundaries);

        for (Gene g : genes.values()) {
            g.sortExonsForEachTranscript();
            String chr = g.getChromosome();
            Boolean key = null;

            if (frStrand != null)
                key = g.getStrand() == '+';

            Map<String, IntervalTree<Gene>> mapForStrand =
                    geneTrees.computeIfAbsent(key, k -> new HashMap<>());
            mapForStrand.computeIfAbsent(chr, c -> new IntervalTree<>()).add(g);
        }
    }

    @FunctionalInterface
    public interface IntervalOp<T> {
        List<T> apply(int start, int end, List<T> out);
    }

    private IntervalTree<Gene> getTree(String chr, Boolean frStrand) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(frStrand);
        if (strandMap == null) return null;
        return strandMap.get(chr);
    }

    public List<Gene> getGenesByIntervalOperation(
            String chr,
            int start,
            int end,
            Boolean frStrand,
            IntervalOp<Gene> op
    ) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(frStrand);
        if (strandMap == null) return Collections.emptyList();
        IntervalTree<Gene> tree = strandMap.get(chr);
        if (tree == null) return Collections.emptyList();
        return op.apply(start, end, new ArrayList<>());
    }

    public List<Gene> getContainingGenes(String chr, int start, int end, Boolean frStrand) {
        IntervalTree<Gene> tree = getTree(chr, frStrand);
        if (tree == null) return Collections.emptyList();
        return getGenesByIntervalOperation(chr, start, end, frStrand, tree::getIntervalsSpanning);
    }

    public List<Gene> getIncludedGene(String chr, int start, int end, Boolean frStrand) {
        IntervalTree<Gene> tree = getTree(chr, frStrand);
        if (tree == null) return Collections.emptyList();
        return getGenesByIntervalOperation(chr, start, end, frStrand, tree::getIntervalsSpannedBy);
    }

    public List<Gene> getRightNeighbor(String chr, int start, int end, Boolean frStrand) {
        IntervalTree<Gene> tree = getTree(chr, frStrand);
        if (tree == null) return Collections.emptyList();
        return getGenesByIntervalOperation(chr, start, end, frStrand, tree::getIntervalsRightNeighbor);
    }

    public List<Gene> getLeftNeighbor(String chr, int start, int end, Boolean frStrand) {
        IntervalTree<Gene> tree = getTree(chr, frStrand);
        if (tree == null) return Collections.emptyList();
        return getGenesByIntervalOperation(chr, start, end, frStrand, tree::getIntervalsLeftNeighbor);
    }

    public List<Gene> getGenes() {
        return genes.values().stream().toList();
    }
}
