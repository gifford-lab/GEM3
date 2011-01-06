package edu.mit.csail.cgs.datasets.proteinprotein;

public class Link {

    public int geneA, geneB;
    public double score;

    public Link(int a, int b, double s) {
        geneA = a;
        geneB = b;
        score = s;
    }

}