package com.github.valentinrexer;

public class ReadPairRegions {
    Region fistRecordRegion;
    Region lastRecordRegion;

    public ReadPairRegions(ReadPair pair) {
        fistRecordRegion = pair.getFirstRecordRegion();
        lastRecordRegion = pair.getLastRecordRegion();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ReadPairRegions that = (ReadPairRegions) o;
        return that.fistRecordRegion.start() == this.fistRecordRegion.start() &&
                that.fistRecordRegion.end() == this.fistRecordRegion.end() &&
                that.lastRecordRegion.start() == this.lastRecordRegion.start() &&
                that.lastRecordRegion.end() == this.lastRecordRegion.end();
    }

    @Override
    public int hashCode() {
        int result = 17;

        result = 31 * result + fistRecordRegion.start();
        result = 31 * result + fistRecordRegion.end();
        result = 31 * result + lastRecordRegion.start();
        result = 31 * result + lastRecordRegion.end();

        return result;
    }
}
