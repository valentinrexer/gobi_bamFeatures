package com.github.valentinrexer;

import java.util.*;

public class PcrIndexMap {
    private final Map<Boolean, Map<Set<Region>, Integer>> pcrIndexMap;

    public PcrIndexMap() {
        pcrIndexMap = new HashMap<>();
    }

    public int getPcrIndex(Set<Region> regionVector, Boolean strand) {
        if (! pcrIndexMap.containsKey(strand))
            pcrIndexMap.put(strand, new HashMap<>());

        var map = pcrIndexMap.get(strand);

        if (!map.containsKey(regionVector)) {
            map.put(regionVector, 1);
            return 0;
        }
        map.put(regionVector, map.get(regionVector) + 1);
        return map.get(regionVector) - 1;
    }

    public int getPcrIndex(List<Region> regionVectorList, Boolean strand) {
        var regionVector = new HashSet<>(regionVectorList);
        if (!pcrIndexMap.containsKey(strand))
            pcrIndexMap.put(strand, new HashMap<>());

        var map = pcrIndexMap.get(strand);

        if (!map.containsKey(regionVector)) {
            map.put(regionVector, 1);
            return 0;
        }
        map.put(regionVector, map.get(regionVector) + 1);
        return map.get(regionVector) - 1;
    }
}
