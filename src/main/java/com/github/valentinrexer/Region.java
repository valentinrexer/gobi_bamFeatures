package com.github.valentinrexer;

import augmentedTree.Interval;

import java.util.Objects;

public record Region(int start, int end) implements Interval {
    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return end;
    }

    public boolean intersects(Region r) {
        return this.start <= r.end() && this.end >= r.start();
    }


    public Region getIntersectingRegion(Region region) {
        if (!intersects(region)) return null;

        return new Region(Math.max(start, region.start()), Math.min(end, region.end()));
    }

    @Override
    public boolean equals(Object o) {
        if (o == null || getClass() != o.getClass()) return false;
        Region region = (Region) o;
        return end == region.end && start == region.start;
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }
}
