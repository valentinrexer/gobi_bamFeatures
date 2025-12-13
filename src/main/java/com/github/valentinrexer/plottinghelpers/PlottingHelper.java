package com.github.valentinrexer.plottinghelpers;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;

import augmentedTree.IntervalTree;
import com.github.valentinrexer.Gene;
import com.github.valentinrexer.Region;
import com.github.valentinrexer.TreeGtf;

public class PlottingHelper {
    private final TreeGtf treeGtf;

    public PlottingHelper(Path pathToGtf) {
        this.treeGtf = new TreeGtf();
        treeGtf.readInGffFile(pathToGtf, null);
    }

    public HashMap<String, Integer> getGeneLengths() {
        var lenMap = new HashMap<String, Integer>();
        for (Gene gene : treeGtf.getGenes()) {
            int mergedLength = 0;

            var mergedTranscriptome = gene.getMergedTranscriptomeForInterval(new Region(Integer.MIN_VALUE, Integer.MAX_VALUE));
            for (Region mergedBlock : mergedTranscriptome)
                mergedLength += mergedBlock.inclusiveLength();

            lenMap.put(gene.getGeneName(), mergedLength);
        }

        return lenMap;
    }

    public static void main(String[] args) {
        Path gtfPath = Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/plotting/complete_bams/Saccharomyces_cerevisiae.R64-1-1.75.gtf");
        Path homoPath = Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/plotting/complete_bams/Homo_sapiens.GRCh37.75.gtf");
        Path outFile = Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/plotting/data/sacc_gene_len.tsv");
        Path outFileHomo = Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/plotting/data/homo_gene_len.tsv");

        PlottingHelper plottingHelper = new PlottingHelper(homoPath);
        var outMap = plottingHelper.getGeneLengths();

        try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFileHomo.toFile())))){
            for (String geneName : outMap.keySet()) {
                bw.write(geneName + "\t" + outMap.get(geneName) + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
