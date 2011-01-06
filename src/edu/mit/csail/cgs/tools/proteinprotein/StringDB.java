package edu.mit.csail.cgs.tools.proteinprotein;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.proteinprotein.*;
import edu.mit.csail.cgs.tools.utils.Args;

public class StringDB {

    private static int species;
    private static edu.mit.csail.cgs.datasets.proteinprotein.StringDB db;

    public static void main(String[] args) throws Exception {
        String speciesname = args[0];
        String cmd = args[1];
        db = new edu.mit.csail.cgs.datasets.proteinprotein.StringDB();
        species = db.getSpeciesID(speciesname);
        if (cmd.equals("geneid")) {
            getGeneID(args[2]);
        } else if (cmd.equals("genelinks")) {
            getGeneLinks(args[2]);
        } else if (cmd.equals("geneactions")) {
            getGeneActions(args[2]);
        }
    }
    private static void getGeneID(String name) throws SQLException, NotFoundException {
        int id = db.getGeneID(name, species);
        System.out.println(name + " is " + id);
    }
    private static void getGeneLinks(String name) throws SQLException, NotFoundException {
        int id = db.getGeneID(name, species);
        List<Link> links = db.getGeneLinks(id);
        for (Link l : links) {            
            String other = db.getGeneName(l.geneB);
            System.out.println(String.format("%s\t%s\t%f",
                                             name, other, l.score));
                                             
        }
    }
    private static void getGeneActions(String name) throws SQLException, NotFoundException {
        int id = db.getGeneID(name, species);
        List<Action> actions = db.getGeneActions(id);
        for (Action a : actions) {            
            String other = db.getGeneName(a.geneB);
            System.out.println(String.format("%s\t%s\t%f\t%s\t%s",
                                             name, other, a.score,a.mode,a.action));
        }
    }
}
