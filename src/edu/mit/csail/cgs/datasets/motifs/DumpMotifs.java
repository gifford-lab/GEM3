package edu.mit.csail.cgs.datasets.motifs;

import java.util.*;
import edu.mit.csail.cgs.tools.utils.Args;


public class DumpMotifs {

    public static void main(String args[]) throws Exception {
        Collection<WeightMatrix> matrices = Args.parseWeightMatrices(args);
        boolean transfac = true;
        for (WeightMatrix m : matrices) {
            if (transfac) {
                printTransfac(m);
            }
        }
    }
    public static void printTransfac(WeightMatrix m) {
        System.out.println("DE\t" + m.toString());
        for (int i = 0; i < m.matrix.length; i++) {
            System.out.println(String.format("%d\t%d\t%d\t%d\t%d\tX",
                                             i, (int)(100 * m.matrix[i]['A']),
                                             (int)(100 * m.matrix[i]['C']),
                                             (int)(100 * m.matrix[i]['G']),
                                             (int)(100 * m.matrix[i]['T'])));
        }
        System.out.println("XX");
    }


}