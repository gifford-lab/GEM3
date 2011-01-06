package edu.mit.csail.cgs.datasets.proteinprotein;

public class Action {

    public int species;
    public String geneA, geneB;
    public String mode, action;
    public double score;

    public Action(int sp, String a, String b, String m, String ac, double s) {
        species = sp;
        geneA = a;
        geneB = b;
        mode = m;
        action = ac;
        score = s;
    }

}
