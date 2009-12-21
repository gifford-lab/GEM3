package edu.mit.csail.cgs.tools.fragdist;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Add a fragment length distribution to the database
 *
 * AddFragDist --fragdist 'name;version' [--description "blah blah blah"] --file data.txt
 *
 * The input file should have one number per line, starting at an offset of zero.
 * The numbers should range from zero to one, inclusive, and represent the shape of
 * the enrichment peak from a single ChIP-Chip binding event
 */

public class AddFragDist {

    public static void main (String args[]) throws Exception {        
        String fname = null;
        String name, version;
        ArrayList<Double> values = new ArrayList<Double>();
        String desc = Args.parseString(args,"description","");
        String nv = Args.parseString(args,"fragdist",null);
        String pieces[] = nv.split(";");
        if (pieces.length != 2) {
            System.err.println("Must supply --fragdist 'name;version'");
            System.exit(1);
        }       
        name = pieces[0];
        version = pieces[1];

        fname = Args.parseString(args,"file",null);
        BufferedReader reader;
        if (fname == null || fname.equals("-")) {
            reader = new BufferedReader(new InputStreamReader(System.in));
        } else {
            reader = new BufferedReader(new FileReader(new File(fname)));
        }
        
        String line;
        while ((line = reader.readLine()) != null) {
            values.add(Double.parseDouble(line));
        }

        double[] array = new double[values.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = values.get(i);
        }
        
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        loader.getFragDist(name,version,desc,array);


        
    }


}