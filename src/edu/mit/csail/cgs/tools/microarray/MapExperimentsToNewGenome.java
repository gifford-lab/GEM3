package edu.mit.csail.cgs.tools.microarray;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;

import java.io.*;
import java.util.*;
import java.sql.*;

/**
 * Maps all experiments with results for one genome to a new genome.
 * Use this with MapProbesToNewGenome; map the probes for a species 
 * to a new genome and then remap the experiments.
 *
 * MapExperimentsToNewGenome --old "$SC;SGDv1" --new "$SC;SGDv2"
 */

public class MapExperimentsToNewGenome {

    public static void main(String args[]) throws Exception {
        String oldsg = Args.parseString(args,"old",null);
        String newsg = Args.parseString(args,"new",null);
        String oldpieces[] = oldsg.split(";");
        String newpieces[] = newsg.split(";");
        Genome oldg = new Genome(oldpieces[0], oldpieces[1]);
        Genome newg = new Genome(newpieces[0], newpieces[1]);
        int oldid = oldg.getDBID();
        int newid = newg.getDBID();
        
        java.sql.Connection c = 
			DatabaseFactory.getConnection("chipchip");
        Statement stmt = c.createStatement();
        PreparedStatement insert = c.prepareStatement("insert into exptToGenome(experiment,genome) values(?,?)");
        
        ResultSet rs = stmt.executeQuery("select experiment from exptToGenome where genome = " + oldid);
        while (rs.next()) {
            int exptid = rs.getInt(1);
            insert.setInt(1,exptid);
            insert.setInt(2,newid);
            insert.execute();
            System.err.println("Mapped " + exptid);
        }
        c.commit();
        c.close();
    }
    
}
