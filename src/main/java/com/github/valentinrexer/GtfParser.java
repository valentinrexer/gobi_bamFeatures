package com.github.valentinrexer;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GtfParser {
    public static GtfData parse(Path gtfPath, Boolean frStrand) throws IOException {

        GtfData data = new GtfData();

        try (BufferedReader reader = Files.newBufferedReader(gtfPath)) {
            String line;

            while ((line = reader.readLine()) != null) {

                // Skip empty or comment lines
                if (line.isEmpty() || line.startsWith("#"))
                    continue;

                String[] fields = line.split("\t");
                if (fields.length < 9)
                    continue;

                String chr = fields[0];
                String featureType = fields[2];
                int start = Integer.parseInt(fields[3]);
                int end = Integer.parseInt(fields[4]);
                char strand = fields[6].charAt(0);
                String attributeString = fields[8];

                // We only care about transcripts + exons
                if (!featureType.equals("transcript") && !featureType.equals("exon"))
                    continue;


                String geneId = getAttribute("gene_id", attributeString);
                String transcriptId = getAttribute("transcript_id", attributeString);

                if (geneId == null || transcriptId == null)
                    continue;

                String geneName = getAttribute("gene_name", attributeString);
                String biotype = getAttribute("gene_biotype", attributeString);

                // Fetch or create gene
                Gene gene = data.getOrCreateGene(geneId, geneName, biotype, strand, chr);

                // Fetch or create transcript
                Transcript transcript = data.getOrCreateTranscript(transcriptId, geneId, strand);

                // Make sure transcript belongs to gene
                gene.addTranscript(transcript);

                // EXON HANDLING
                if (featureType.equals("exon")) {
                    Exon exon = new Exon(start, end, transcriptId);
                    transcript.addExon(exon);
                }
            }
        }

        data.buildIntervalTrees(frStrand != null);

        return data;
    }

    private static String getAttribute(String attribute, String s) {
        Pattern pattern = Pattern.compile(attribute + "\\s+\"([^\"]+)\"");
        Matcher matcher = pattern.matcher(s);
        if (matcher.find()) {
            return matcher.group(1);
        } else {
            System.out.println("Wasn't able to extract '" + attribute + "' for line object " + s);
            return null;
        }
    }

    private static Map<String, String> parseAttributes(String attr) {
        Map<String, String> map = new HashMap<>();

        // attributes are key "value"; key "value";
        String[] tokens = attr.split(";");
        for (String tok : tokens) {
            tok = tok.trim();
            if (tok.isEmpty())
                continue;

            int spaceIdx = tok.indexOf(' ');
            if (spaceIdx < 0)
                continue;

            String key = tok.substring(0, spaceIdx).trim();
            String val = tok.substring(spaceIdx + 1).trim();

            // Strip quotes
            val = val.replace("\"", "").trim();

            if (!key.isEmpty() && !val.isEmpty())
                map.put(key, val);
        }

        return map;
    }
}
