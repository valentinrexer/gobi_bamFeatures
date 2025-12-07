package com.github.valentinrexer;

public class Runner {
    public static void main(String[] args) {
        Region r1 = new Region(10, 20);
        Region r2 = new Region(9, 25);

        System.out.println(r1.getIntersectingRegion(r2));
    }
}
