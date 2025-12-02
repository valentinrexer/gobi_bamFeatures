package com.github.valentinrexer;

import htsjdk.samtools.*;

import java.nio.file.Path;
import java.nio.file.Paths;
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
}
