package com.github.valentinrexer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader;

import java.nio.file.Path;
import java.util.Iterator;

public class SamFile {
    private final SamReader reader;

    public SamFile(Path samFile) {
        reader = SamReaderFactory
                .makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(samFile.toFile());
    }

    public SamReader getReader() {
        return reader;
    }

    public Iterator<SAMRecord> iterator() {
        return reader.iterator();
    }

    public SAMFileHeader.SortOrder sortOrder() {
        SAMFileHeader header = reader.getFileHeader();
        return header.getSortOrder();
    }

}
