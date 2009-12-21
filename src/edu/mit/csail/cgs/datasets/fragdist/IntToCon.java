package edu.mit.csail.cgs.datasets.fragdist;

import java.util.*;
import java.io.*;

/**
 * IntToCon process Agilent Bioanalyzer data by
 * turning intensity measurements into concentrations.
 * Since dye incorporates into the length of a molecule,
 * the number of molecules times their length gives intensity.
 *
 * IntToCon reads from stdin and writes to stdout.  The
 * intput lines are tab separated lengths and intensities.
 * 
 */

public class IntToCon {

    public static void main(String args[]) throws IOException {
        String line;
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        while ((line = reader.readLine()) != null) {
            if (line.matches("[\\.\\d]+\\t[\\d\\.]+")) {
                String pieces[] = line.split("\\t");                
                Double len = Double.parseDouble(pieces[0]);
                Double intensity = Double.parseDouble(pieces[1]);
                System.out.println(pieces[0] + "\t" + (intensity / len));
            }
        }
    }

}
