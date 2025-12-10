package com.github.valentinrexer;

import augmentedTree.IntervalTree;
import java.util.*;

@Deprecated
public class GtfData {

    private final Map<String, Gene> genesById = new HashMap<>();
    private final Map<String, Transcript> transcriptsById = new HashMap<>();
    private final Map<Boolean, Map<String, IntervalTree<Gene>>> geneTrees = new HashMap<>();

    public Gene getOrCreateGene(String geneId, String geneName, String biotype, char strand, String chr) {
        return genesById.computeIfAbsent(geneId, id -> new Gene(id, geneName, biotype, strand, chr));
    }

    public Transcript getOrCreateTranscript(String txId, String geneId, char strand) {
        return transcriptsById.computeIfAbsent(txId, id -> new Transcript(txId, geneId, strand));
    }

    public void buildIntervalTrees(boolean strandedExperiment) {
        genesById.values().forEach(Gene::computeBoundaries);

        for (Gene g : genesById.values()) {
            g.sortExonsForEachTranscript();
            String chr = g.getChromosome();
            Boolean key = strandedExperiment ? strandKey(g.getStrand()) : null;
            Map<String, IntervalTree<Gene>> mapForStrand =
                    geneTrees.computeIfAbsent(key, k -> new HashMap<>());
            mapForStrand.computeIfAbsent(chr, c -> new IntervalTree<>()).add(g);
        }
    }

    private Boolean strandKey(char strand) {
        if (strand == '+') return Boolean.TRUE;
        if (strand == '-') return Boolean.FALSE;
        return null;
    }

    public List<Gene> getContainingGenes(String chr, int start, int end, Boolean sense) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(sense);
        if (strandMap == null) return Collections.emptyList();
        IntervalTree<Gene> tree = strandMap.get(chr);
        if (tree == null) return Collections.emptyList();
        return tree.getIntervalsSpanning(start, end, new ArrayList<>());
    }

    public List<Gene> getIncludedGene(String chr, int start, int end, Boolean sense) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(sense);
        if (strandMap == null) return Collections.emptyList();
        IntervalTree<Gene> tree = strandMap.get(chr);
        if (tree == null) return Collections.emptyList();
        return tree.getIntervalsSpannedBy(start, end, new ArrayList<>());
    }

    public List<Gene> getRightNeighbor(String chr, int start, int end, Boolean sense) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(sense);
        if (strandMap == null) return Collections.emptyList();
        IntervalTree<Gene> tree = strandMap.get(chr);
        if (tree == null) return Collections.emptyList();
        return tree.getIntervalsRightNeighbor(start, end, new ArrayList<>());
    }

    public List<Gene> getLeftNeighbor(String chr, int start, int end, Boolean sense) {
        Map<String, IntervalTree<Gene>> strandMap = geneTrees.get(sense);
        if (strandMap == null) return Collections.emptyList();
        IntervalTree<Gene> tree = strandMap.get(chr);
        if (tree == null) return Collections.emptyList();
        return tree.getIntervalsLeftNeighbor(start, end, new ArrayList<>());
    }
}
