/*
 * Created on Sep 14, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.chipchip;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

/**
 * Usage:
 * java edu.mit.csail.cgs.tools.chipchip.FindUniqueProbeInformation SGDv1 "Sc Gcn4:Sc:YPD vs WCE:Sc:YPD" "median linefit" 
 *
 * Prints information about each probe in the specified experiment as mapped to the specified Genome.
 * Lines starting with a non-whitespace character give probe information (probe id, name, expt id, cy5 value, cy3 value).
 * Lines starting with a tab give the probe location and score for the preceeding probe.
 */

public class FindUniqueProbeInformation {
    
    public static void main(String[] args) { 
        try {
            Genome g = Organism.findGenome(args[0]);
            String exptName = args[1];
            String exptVersion = args[2];
            ChipChipLocator loc = new ChipChipLocator(g, exptName, exptVersion);
            
            FindUniqueProbeInformation info = new FindUniqueProbeInformation(loc, g);
            info.countUniqueProbeIDs();
            
            Vector<Integer> dbids = info.findExperimentIDs();
            for(int exptID : dbids) { 
                info.outputExperimentData(exptID);
            }
            
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    private Genome genome;
    private ChipChipLocator locator;
    private Map<Integer,Vector<StrandedRegion>> probeLocs;
    private Map<Integer,Vector<Integer>> probeCounts, probeBitScores;
    private java.sql.Connection cc_cxn;

    public FindUniqueProbeInformation(ChipChipLocator loc, Genome g) throws SQLException { 
        locator = loc;
        genome = g;
        cc_cxn = DatabaseFactory.getConnection("chipchip");
    }
    
    public Vector<Integer> findExperimentIDs() throws SQLException {
        Vector<Integer> dbids = new Vector<Integer>();
        Statement s = cc_cxn.createStatement();
        String query = String.format(
                "select id from experiment e, exptToGenome e2g where e.name='%s' and e.version='%s' " +
                "and e.id=e2g.experiment and e2g.genome=%d", locator.name,
                locator.version, genome.getDBID());
        
        ResultSet rs = s.executeQuery(query);
        while(rs.next()) { 
            dbids.add(rs.getInt(1));
        }
        rs.close();
        
        s.close();
        
        return dbids;
    }
    
    public void outputExperimentData(int exptID) throws SQLException {
        System.out.println("Printing Experimental Data...");
        Statement s = cc_cxn.createStatement();
        
        String query = String.format(
                "select d.channelone, d.channeltwo, pd.probename, pd.probeid, pd.id from data d, probedesign pd " +
                "where d.probe=pd.id and d.experiment=%d",
                exptID);

        ResultSet rs = s.executeQuery(query);
        while(rs.next()) { 
            double ip = rs.getDouble(1);
            double wce = rs.getDouble(2);
            String pname = rs.getString(3);
            String pid = rs.getString(4);
            int dbid = rs.getInt(5);

            if(probeLocs.containsKey(dbid)) { 
                System.out.println(String.format(
                        "%s\t%s\t%d\t%f\t%f", pid, pname, exptID, ip, wce));

                for(int i = 0; i < probeLocs.get(dbid).size(); i++) { 
                    String locString = probeLocs.get(dbid).get(i).getLocationString();
                    char strand = probeLocs.get(dbid).get(i).getStrand();
                    int counts = probeCounts.get(dbid).get(i);
                    int bits = probeBitScores.get(dbid).get(i);
                    System.out.println(String.format("\t%s\t%c\t%d\t%d", locString, strand, counts, bits));
                }
            } else { 
                System.out.println(String.format(
                        "%s\t%s\t%d\t%f\t%f\t***", pid, pname, exptID, ip, wce));                
            }
        }
        rs.close();
        
        s.close();
    }
    
    public void countUniqueProbeIDs() throws SQLException {
        System.out.println("Counting Unique Probes...");
        java.sql.Connection cc_cxn = DatabaseFactory.getConnection("chipchip");
        Statement s = cc_cxn.createStatement();
        
        String query = String.format(
                "select pl.id, c.name, pl.startpos, pl.stoppos, pl.strand, pl.loccount, pl.bitscore " +
                "from probelocation pl, core.chromosome c " +
                "where pl.chromosome=c.id and c.genome=%d",
                genome.getDBID());

        probeLocs = new HashMap<Integer,Vector<StrandedRegion>>();
        probeCounts = new HashMap<Integer,Vector<Integer>>();
        probeBitScores = new HashMap<Integer,Vector<Integer>>();
        
        ResultSet rs = s.executeQuery(query);
        while(rs.next()) { 
            int dbid = rs.getInt(1);
            if(!probeLocs.containsKey(dbid)) { 
                probeLocs.put(dbid, new Vector<StrandedRegion>());
                probeCounts.put(dbid, new Vector<Integer>());
                probeBitScores.put(dbid, new Vector<Integer>());
            }
            
            String chrom = rs.getString(2);
            int start = rs.getInt(3);
            int end = rs.getInt(4);
            char strand = rs.getString(5).charAt(0);
            int count = rs.getInt(6);
            int bits = rs.getInt(7);
            
            StrandedRegion strreg = new StrandedRegion(genome, chrom, start, end, strand);
            probeLocs.get(dbid).add(strreg);
            probeCounts.get(dbid).add(count);
            probeBitScores.get(dbid).add(bits);
        }
        rs.close();
        s.close();
        DatabaseFactory.freeConnection(cc_cxn);
        
        System.out.println(String.format("Counted %d probes.", probeLocs.size()));
    }
}
