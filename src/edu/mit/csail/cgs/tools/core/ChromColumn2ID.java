package edu.mit.csail.cgs.tools.core;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import java.io.*;

/**
 * Replaces a chromosome name with its ID in a column of tab-delimited text.
 *
 * eg,
 * cat PPG_D5_WCE-iHoxc6_1.bowtie.align.gz | gunzip | java edu.mit.csail.cgs.projects.readdb.BowtieToReadDB | \
 *   java edu.mit.csail.cgs.tools.core.ChromColumn2ID --species "$MM;mm8" --column 0
 *
 */

public class ChromColumn2ID {
    
    public static void main(String args[]) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        int field = Args.parseInteger(args,"column",-1);
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line;
        StringBuffer out = new StringBuffer();
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");
            pieces[field] = Integer.toString(genome.getChromID(pieces[field].replaceAll("^chr","")));
            out.delete(0,out.length());
            out.append(pieces[0]);
            for (int i = 1; i < pieces.length ; i++) {
                out.append("\t" + pieces[i]);
            }

            System.out.println(out.toString());
        }
    }
}