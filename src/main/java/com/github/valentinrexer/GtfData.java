package com.github.valentinrexer;

import augmentedTree.IntervalTree;
import java.util.*;

public class GtfData {

    private final Map<String, Gene> genesById = new HashMap<>();
    private final Map<String, Transcript> transcriptsById = new HashMap<>();

    private final Map<String, IntervalTree<Gene>> geneSenseTrees = new HashMap<>();
    private final Map<String, IntervalTree<Gene>> geneAntiSenseTrees = new HashMap<>();

    public Gene getOrCreateGene(String geneId, String geneName, String biotype, char strand, String chr) {
        return genesById.computeIfAbsent(
                geneId, id -> new Gene(id, geneName, biotype, strand, chr)
        );
    }

    public Transcript getOrCreateTranscript(String txId, String geneId, char strand) {
        return transcriptsById.computeIfAbsent(
                txId, id -> new Transcript(txId, geneId, strand)
        );
    }

    public void buildIntervalTrees() {
        genesById.values().forEach(Gene::computeBoundaries);

        for (Gene g : genesById.values()) {
            g.sortExonsForEachTranscript();

            String chr = g.getChromosome();

            if (g.getStrand() == '+') {
                geneSenseTrees
                        .computeIfAbsent(chr, c -> new IntervalTree<>())
                        .add(g);
            } else {
                geneAntiSenseTrees
                        .computeIfAbsent(chr, c -> new IntervalTree<>())
                        .add(g);
            }
        }
    }

    public List<Gene> getOverlappingGenesSense(String chr, int start, int end) {
        IntervalTree<Gene> tree = geneSenseTrees.get(chr);
        if (tree == null) return Collections.emptyList();

        return tree.getIntervalsIntersecting(start, end, new ArrayList<>());
    }

    public List<Gene> getOverlappingGenesAntisense(String chr, int start, int end) {
        IntervalTree<Gene> tree = geneAntiSenseTrees.get(chr);
        if (tree == null) return Collections.emptyList();

        return tree.getIntervalsIntersecting(start, end, new ArrayList<>());
    }

    public List<Gene> getRightNeighborSense(String chr, int pos) {
        IntervalTree<Gene> tree = geneSenseTrees.get(chr);
        if (tree == null) return Collections.emptyList();

        return tree.getIntervalsRightNeighbor(pos, pos, new ArrayList<>());
    }

    public List<Gene> getLeftNeighborSense(String chr, int pos) {
        IntervalTree<Gene> tree = geneSenseTrees.get(chr);
        if (tree == null) return Collections.emptyList();

        return tree.getIntervalsLeftNeighbor(pos, pos, new ArrayList<>());
    }
}
