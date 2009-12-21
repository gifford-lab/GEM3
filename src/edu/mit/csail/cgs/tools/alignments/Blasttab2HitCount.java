package edu.mit.csail.cgs.tools.alignments;

import java.util.*;
import java.io.*;

/* command line program that takes blasttab input on stdin.
   Produces on stdout lines to load into a mysql table with schema:
+------------+--------------+------+-----+---------+-------+
| Field      | Type         | Null | Key | Default | Extra |
+------------+--------------+------+-----+---------+-------+
| chrom      | varchar(255) | YES  | MUL | NULL    |       | 
| chromStart | int(10)      | YES  |     | NULL    |       | 
| chromEnd   | int(10)      | YES  |     | NULL    |       | 
| score      | double       | YES  |     | NULL    |       | 
+------------+--------------+------+-----+---------+-------+

 that shows how many hits against the target genome there are for
 each region in the source genome.  The chromosomes and coordinates in the
 output are from the source genome and the score is the number of hits
 in the target genome

 blasttab format is 
  0) query id
 1) target id
 2) match percent
 3) match len
 4)
 5) 
 6) qstart
 7) qend
 8) tstart
 9) tend
10) e-value
11) score

2087597 chr14   100.00  60      0       0       1       60      37301   37360   6e-28    119
2087597 chr14   100.00  13      0       0       44      56      306606  306594  6.7     26.3
2087597 chr14   100.00  13      0       0       13      25      576177  576165  6.7     26.3
2087597 chr14   100.00  13      0       0       42      54      666356  666344  6.7     26.3

*/

public class Blasttab2HitCount {

    public static void main(String args[]) throws IOException {
        Map<String,List<Integer>> map = new HashMap<String,List<Integer>>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = reader.readLine()) != null) {
            String[] pieces = line.split("\\t");
            if (!map.containsKey(pieces[0])) {
                map.put(pieces[0], new ArrayList<Integer>());
            }
            List<Integer> list = map.get(pieces[0]);
            for (int pos = Integer.parseInt(pieces[6]); pos <= Integer.parseInt(pieces[7]); pos++) {
                if (pos < list.size()) {
                    list.set(pos, list.get(pos) + 1);
                } else {
                    for (int j = list.size(); j < pos; j++) {
                        list.add(0);
                    }
                    list.add(1);
                }
            }
        }
        for (String chrom : map.keySet()) {
            List<Integer> list = map.get(chrom);
            for (int pos = 0; pos < list.size(); pos++) {
                int count = list.get(pos);
                int j = pos;
                while (j + 1 < list.size() && list.get(j+1) == count) {
                    j++;
                }
                if (count != 0) {
                    System.out.println(String.format("%s\t%d\t%d\t%d",
                                                     chrom,
                                                     pos,
                                                     j,
                                                     count));
                }
                pos = j;
            }

        }

    }

}