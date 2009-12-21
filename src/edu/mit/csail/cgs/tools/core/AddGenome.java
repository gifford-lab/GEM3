package edu.mit.csail.cgs.tools.core;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * AddGenome --species "mus musculus" --genome "mm22"
 */

public class AddGenome {

    public static void main(String args[]) throws Exception {
        String species = Args.parseString(args,"species",null);
        String genome = Args.parseString(args,"genome",null);
        if (species == null) {
            System.err.println("Must supply --species");
            System.exit(1);
        }
        if (genome == null) {
            System.err.println("Must supply --genome");
            System.exit(1);
        }
        try {
            Organism o = new Organism(species);
            try {
                Genome g = o.getGenome(genome);
                System.err.println("Genome " + genome + " already exists.");
                System.exit(2);
            } catch (NotFoundException e) {
                o.insertGenome(genome);
            }
        } catch (NotFoundException e) {
            System.err.println("Species " + species + " doesn't exist.  Please create with AddSpecies");
            System.exit(2);
        }


    }

}