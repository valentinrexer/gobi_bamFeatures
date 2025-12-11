package com.github.valentinrexer.utils;

import com.github.valentinrexer.Region;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;

public class BamFeatureUtils {
    public static List<Region> mergeVector(List<Region> vector) {
        if (vector.isEmpty()) return List.of();

        // remove duplicates
        vector = new ArrayList<>(new LinkedHashSet<>(vector));

        // sort
        vector.sort(Comparator.comparingInt(Region::start)
                .thenComparingInt(Region::end));

        var merged = new ArrayList<Region>();
        Region prev = vector.get(0);

        for (int i = 1; i < vector.size(); i++) {
            Region next = vector.get(i);

            if (prev.end() >= next.start() - 1) {
                prev = new Region(prev.start(), Math.max(prev.end(), next.end()));
            } else {
                merged.add(prev);
                prev = next;
            }
        }
        merged.add(prev);

        return merged;
    }

}
