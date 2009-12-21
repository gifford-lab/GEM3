package edu.mit.csail.cgs.tools.core;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * AddSpecies --species "mus musculus"
 */

public class AddSpecies {

    public static void main(String args[]) throws Exception {
        String species = Args.parseString(args,"species",null);
        if (species == null) {
            System.err.println("Must supply --species");
            System.exit(1);
        }
        try {
            Organism o = new Organism(species);
            System.err.println("Species " + species + " already exists");
            System.exit(2);
        } catch (NotFoundException e) {
            Organism.insertOrganism(species);
        }
    }

}