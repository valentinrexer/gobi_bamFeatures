package com.github.valentinrexer;

public class Exon {
    private final int start;
    private final int end;
    private final String transcriptId;

    public Exon(int start, int end, String transcriptId) {
        this.start = start;
        this.end = end;
        this.transcriptId = transcriptId;
    }

    public int getStart() { return start; }
    public int getEnd() { return end; }
    public String getTranscriptId() { return transcriptId; }

    @Override
    public String toString() {
        return "Exon{" + start + "-" + end + ", tx=" + transcriptId + "}";
    }
}

