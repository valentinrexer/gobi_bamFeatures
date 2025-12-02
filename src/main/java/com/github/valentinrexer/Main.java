package com.github.valentinrexer;

import htsjdk.samtools.*;

import java.nio.file.Paths;

public class Main {
    public static void main(String[] args) {
        SamFile sam = new SamFile(Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/BamFeatures/complete_bams/nookaew_cm.bam"));

        for (SAMRecord record : sam.getReader())
            System.out.println(record);
    }
}