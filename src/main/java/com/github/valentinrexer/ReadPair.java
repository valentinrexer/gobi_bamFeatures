package com.github.valentinrexer;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class ReadPair {
    private final SAMRecord firstRecord;
    private final SAMRecord lastRecord;
    private final List<Region> regionVector;

    public ReadPair(SAMRecord firstRecord, SAMRecord lastRecord) {
        this.firstRecord = firstRecord;
        this.lastRecord = lastRecord;

        List<Region> regionVector = getRegionVector(firstRecord);
        regionVector.addAll(getRegionVector(lastRecord));
        this.regionVector = mergeRegions(regionVector);
    }

    public void process(GtfData gtfData) {
        int mm = getMismatches();
        int clipped = getTotalClipped();
        int nSplit = getNSplit(firstRecord) + getNSplit(lastRecord);

        System.out.println(firstRecord.getReadName() + " " + getCandidateGenes(gtfData).size());

    }

    private List<Gene> getCandidateGenes(GtfData gtfData) {
        String chr = firstRecord.getReferenceName();
        Region region = getRVBoundaries();
        return gtfData.getOverlappingGenesSense(chr, region.start(), region.end());
    }

    private Region getRVBoundaries() {
        int minStart = Integer.MAX_VALUE;
        int maxEnd = Integer.MIN_VALUE;

        for (Region region : regionVector) {
            minStart = Math.min(minStart, region.start());
            maxEnd = Math.max(maxEnd, region.end());
        }

        return new Region(minStart, maxEnd, false);
    }

    private static List<Region> mergeRegions(List<Region> regions) {
        if (regions.isEmpty()) return regions;

        regions.sort(Comparator
                .comparingInt(Region::start)
                .thenComparingInt(Region::end));

        List<Region> merged = new ArrayList<>();
        Region prev = regions.getFirst();

        for (int i = 1; i < regions.size(); i++) {
            Region curr = regions.get(i);

            if (prev.isIntronic() == curr.isIntronic() &&
                    curr.start() <= prev.end() + 1) {

                prev = new Region(
                        prev.start(),
                        Math.max(prev.end(), curr.end()),
                        prev.isIntronic()
                );

            } else {
                merged.add(prev);
                prev = curr;
            }
        }

        merged.add(prev);
        return merged;
    }

    private static List<Region> getRegionVector(SAMRecord record) {
        List<Region> regions = new ArrayList<>();

        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        var pointer = record.getAlignmentStart();

        for (CigarElement cigarElement : cigarElements) {
            var length = cigarElement.getLength();
            CigarOperator operator = cigarElement.getOperator();

            switch (operator) {
                case M:
                case EQ:
                case X:
                case D:
                    regions.add(new Region(pointer, pointer + length - 1, false));
                    pointer += length;
                    break;

                case N:
                    regions.add(new Region(pointer, pointer + length - 1, true));
                    pointer += length;
                    break;

                case I:
                case S:
                case H:
                case P:
                    break;

                default:
                    throw new IllegalStateException("Unknown operator: " + operator);
            }
        }

        return regions;
    }

    private static int getNSplit(SAMRecord record) {
        int nCount = 0;

        for (CigarElement element : record.getCigar().getCigarElements()) {
            if (element.getOperator() == CigarOperator.N) nCount++;
        }

        return nCount;
    }

    private int getMismatches() {
        Integer nm1 = firstRecord.getIntegerAttribute("NM");
        Integer nm2 = lastRecord.getIntegerAttribute("NM");

        if (nm1 == null) nm1 = 0;
        if (nm2 == null) nm2 = 0;

        return nm1 + nm2;
    }

    private int getTotalClipped() {
        return getClippingCount(firstRecord) + getClippingCount(lastRecord);
    }

    private static int getClippingCount(SAMRecord record) {
        List<CigarElement> cigar = record.getCigar().getCigarElements();
        int clipped = 0;

        for (CigarElement element : cigar) {
            CigarOperator operator = element.getOperator();

            if (operator == CigarOperator.S || operator == CigarOperator.H) {
                clipped += element.getLength();
            }
        }

        return clipped;
    }

    @Override
    public String toString() {
        return firstRecord.getReadName() + ", " + lastRecord.getReadName();
    }
}
