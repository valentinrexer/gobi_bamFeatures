package com.github.valentinrexer;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

import javax.swing.*;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Objects;

public class Main {
    public static void main(String[] args) throws IOException {

        // One pending-record map per chromosome
        Map<String, Map<String, SAMRecord>> pendingPerChromosome = new HashMap<>();

        GtfData gtfData = GtfParser.parse(
                Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/BamFeatures/Saccharomyces_cerevisiae.R64-1-1.75.gtf"),
        null);
        SamFile sam = new SamFile(
                Paths.get("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/bam_data/BamFeatures/y.ns.2.bam")
        );

        var c = 0;

        for (Iterator<SAMRecord> it = sam.iterator(); it.hasNext(); ) {
            SAMRecord record = it.next();

            if (!record.getReadPairedFlag()) continue;
            if (record.getReadUnmappedFlag()) continue;
            if (record.getMateUnmappedFlag()) continue;
            if (record.isSecondaryOrSupplementary()) continue;
            if (!Objects.equals(record.getReferenceIndex(), record.getMateReferenceIndex())) continue;
            if (record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag()) continue;

            String chr = record.getReferenceName();
            pendingPerChromosome.putIfAbsent(chr, new HashMap<>());
            Map<String, SAMRecord> pendingRecords = pendingPerChromosome.get(chr);

            String readName = record.getReadName();

            if (!pendingRecords.containsKey(readName)) {
                pendingRecords.put(readName, record);
                continue;
            }

            SAMRecord pendingRecord = pendingRecords.remove(readName);
            SAMRecord first, second;

            if (record.getFirstOfPairFlag() && pendingRecord.getSecondOfPairFlag()) {
                first = record;
                second = pendingRecord;
            } else if (record.getSecondOfPairFlag() && pendingRecord.getFirstOfPairFlag()) {
                first = pendingRecord;
                second = record;
            } else {
                continue;
            }

            ReadPair pair = new ReadPair(first, second);
            pair.process(gtfData, null);
        }
        System.out.println(c);
    }
}
