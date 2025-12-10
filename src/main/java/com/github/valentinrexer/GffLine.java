package com.github.valentinrexer;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GffLine {
    /*
    Holds all Information on exactly one Line in the Gff File
     */

    private final String seqId;
    private final String source;
    private final String type;
    private final int start;
    private final int end;
    private final String score;
    private final char strand;
    private final String phase;
    private final String attributes;

    public GffLine(String line) {
        String[] parts = line.split("\t");
        if (parts.length != 9) {
            throw new IllegalArgumentException("Ungültige GFF-Zeile: " + line);
        }
        this.seqId = parts[0];
        this.source = parts[1];
        this.type = parts[2];
        this.start = Integer.parseInt(parts[3]);
        this.end = Integer.parseInt(parts[4]);
        this.score = parts[5];
        this.strand = parts[6].charAt(0);
        this.phase = parts[7];
        this.attributes = parts[8];
    }

    public GffLine(String seqId, String source, String type, int start, int end,
                   String score, char strand, String phase, String attributes) {
        this.seqId = seqId;
        this.source = source;
        this.type = type;
        this.start = start;
        this.end = end;
        this.score = score;
        this.strand = strand;
        this.phase = phase;
        this.attributes = attributes;
    }

    public String getSeqId() { return seqId; }
    public String getSource() { return source; }
    public String getType() { return type; }
    public int getStart() { return start; }
    public int getEnd() { return end; }
    public String getScore() { return score; }
    public char getStrand() { return strand; }
    public String getPhase() { return phase; }
    public String getAttributes() { return attributes; }

    public String getAttribute(String attribute) {
        Pattern pattern = Pattern.compile(attribute + "\\s+\"([^\"]+)\"");
        Matcher matcher = pattern.matcher(this.attributes);
        if (matcher.find()) {
            return matcher.group(1);
        } else {
            return null;
        }
    }

    @Override
    public String toString() {
        return String.join("\t",
                seqId,
                source,
                type,
                String.valueOf(start),
                String.valueOf(end),
                score,
                String.valueOf(strand),
                phase,
                attributes);
    }
}
