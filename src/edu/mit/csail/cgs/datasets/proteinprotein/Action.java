package edu.mit.csail.cgs.datasets.proteinprotein;

public class Action {

    public int geneA, geneB;
    public String mode, action;
    public double score;

    public Action(int a, int b, String m, String ac, double s) {
        geneA = a;
        geneB = b;
        mode = m;
        action = ac;
        score = s;
    }

}
