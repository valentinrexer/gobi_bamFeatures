package com.github.valentinrexer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.apache.commons.cli.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

public class Main {
    public static void main(String[] args) throws IOException {
        Options options = new Options();

        options.addOption(Option.builder("gtf")
                .hasArg()
                .argName("gtf_file")
                .required(true)
                .desc("GTF annotation file")
                .build());

        options.addOption(Option.builder("bam")
                .hasArg()
                .argName("bam_file")
                .required(true)
                .desc("Input BAM file")
                .build());

        options.addOption(Option.builder("o")
                .longOpt("output")
                .hasArg()
                .argName("output_tsv")
                .required(true)
                .desc("Output TSV file")
                .build());

        options.addOption(Option.builder("frstrand")
                .hasArg()
                .argName("true/false")
                .required(false)
                .desc("FR-stranded flag")
                .build());

        CommandLineParser cliParser = new DefaultParser();
        CommandLine cmd;

        try {
            cmd = cliParser.parse(options, args);
        } catch (ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("bamfeatures", options, true);
            System.err.println("Error: " + e.getMessage());
            return;
        }

        Logger logger = Logger.getLogger(Main.class.getName());

        Path gtfPath = Paths.get(cmd.getOptionValue("gtf"));
        Path bamPath = Paths.get(cmd.getOptionValue("bam"));
        Path outPath = Paths.get(cmd.getOptionValue("o"));

        Boolean frStrand = null;
        if (cmd.hasOption("frstrand")) {
            frStrand = Boolean.parseBoolean(cmd.getOptionValue("frstrand"));
        }

        TreeGtf treeGtf = new  TreeGtf();
        treeGtf.readInGffFile(gtfPath, frStrand);

        SamReader sam = SamReaderFactory
                .makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(bamPath.toFile());

        String currentChromosome = "";
        Map<String, SAMRecord> pendingRecords = new HashMap<>();
        PcrIndexMap pcrIndexMap = new PcrIndexMap();

        try (BufferedWriter writer = Files.newBufferedWriter(outPath)) {
            for (SAMRecord record : sam) {
                if (!record.getReadPairedFlag()) continue;
                if (record.getReadUnmappedFlag()) continue;
                if (record.getMateUnmappedFlag()) continue;
                if (record.isSecondaryOrSupplementary()) continue;
                if (!record.getReferenceName().equals(record.getMateReferenceName())) continue;
                if (record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag()) continue;

                String chr = record.getReferenceName();
                if (!chr.equals(currentChromosome)) {
                    pendingRecords.clear();
                    currentChromosome = chr;
                }

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
                String result = pair.process(treeGtf, frStrand, pcrIndexMap);

                writer.write(result);
                writer.newLine();
                logger.info("Processed " + record.getReadName());
            }
        }
    }
}
