package edu.mit.csail.cgs.datasets.proteinprotein;

public class Link {

    public int species;
    public String geneA, geneB;
    public double score;

    public Link(int sp, String a, String b, double s) {
        species = sp;
        geneA = a;
        geneB = b;
        score = s;
    }

}