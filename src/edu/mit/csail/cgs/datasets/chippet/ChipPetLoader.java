/*
 * Created on May 16, 2007
 */
package edu.mit.csail.cgs.datasets.chippet;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class ChipPetLoader implements edu.mit.csail.cgs.utils.Closeable {
    
    public static void main(String[] args) { 
        try {
            ChipPetLoader loader = new ChipPetLoader();
            Collection<ChipPetExpt> expts = loader.loadAllExperiments();
            for(ChipPetExpt expt : expts) {
                Collection<Genome> gens = loader.loadExperimentGenomes(expt);
                System.out.print("(" + expt.getDBID() + ") " + expt.getName() + " --> ");
                for(Genome g : gens) { 
                    System.out.print(g.getVersion() + " ");
                }
                System.out.println();
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
	
    private java.sql.Connection cxn;
    private PreparedStatement loadByRegion;
    
    public ChipPetLoader() throws SQLException { 
        cxn = DatabaseFactory.getConnection("chippet");
        loadByRegion = cxn.prepareStatement("select startpos, stoppos, strand, peakOverlap from " +
                "chippetdata where expt=? and chromosome=? and ((startpos <= ? and stoppos >= ?) or " +
                "(startpos >= ? and startpos <= ?))");
    }
    
    public Collection<Genome> loadExperimentGenomes(ChipPetExpt expt) throws SQLException { 
        LinkedList<Genome> genomes = new LinkedList<Genome>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select genome from chippetToGenome where expt=" + expt.getDBID());
        while(rs.next()) { 
            int gid = rs.getInt(1);
            try {
                Genome g = Organism.findGenome(gid);
                genomes.add(g);
            } catch (NotFoundException e) {
                e.printStackTrace();
            }
        }
        rs.close();
        s.close();
        return genomes;
    }
    
    public Collection<ChipPetExpt> loadAllExperiments() throws SQLException { 
        PreparedStatement ps = ChipPetExpt.createLoadAll(cxn);
        LinkedList<ChipPetExpt> expts = new LinkedList<ChipPetExpt>();
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            expts.addLast(new ChipPetExpt(rs));
        }
        rs.close();
        ps.close();
        
        return expts;
    }
    
    public ChipPetExpt loadExperiment(String name) throws SQLException { 
        PreparedStatement ps = ChipPetExpt.createLoadByName(cxn);
        ps.setString(1, name);
        ResultSet rs = ps.executeQuery();
        ChipPetExpt expt = null;
        if(rs.next()) { 
            expt = new ChipPetExpt(rs); 
        }
        rs.close();
        ps.close();
        
        if(expt == null) { throw new IllegalArgumentException(name); }
        return expt;
    }
    
    public ChipPetExpt loadExperiment(int dbid) throws SQLException { 
        PreparedStatement ps = ChipPetExpt.createLoadByDBID(cxn);
        ps.setInt(1, dbid);
        ResultSet rs = ps.executeQuery();
        ChipPetExpt expt = null;
        if(rs.next()) { 
            expt = new ChipPetExpt(rs); 
        }
        rs.close();
        ps.close();
        
        if(expt == null) {
        	String err = String.format("No such ChipPet Experiment %d", dbid);
        	throw new IllegalArgumentException(err); 
        }
        return expt;
    }
    
    public Vector<ChipPetDatum> loadAllData(Genome g, ChipPetExpt expt) throws SQLException { 
        Vector<ChipPetDatum> data = new Vector<ChipPetDatum>();
        PreparedStatement ps = ChipPetDatum.prepareLoadByExptAndGenome(cxn);
        ps.setInt(1, expt.getDBID());
        ps.setInt(2, g.getDBID());
        
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            ChipPetDatum datum = new ChipPetDatum(rs, g, expt);
            data.add(datum);
        }
        rs.close();
        ps.close();
        return data;
    }
    
    public Collection<ChipPetDatum> loadIntervals(ChipPetExpt expt, Region r) throws SQLException {
        Genome g = r.getGenome();
        LinkedList<ChipPetDatum> intvs = new LinkedList<ChipPetDatum>();
        synchronized (loadByRegion) {
            loadByRegion.setInt(1, expt.getDBID());
            loadByRegion.setInt(2, g.getChromID(r.getChrom()));
            
            loadByRegion.setInt(3, r.getStart());
            loadByRegion.setInt(4, r.getStart());
            loadByRegion.setInt(5, r.getStart());
            loadByRegion.setInt(6, r.getEnd());
        
            ResultSet rs = loadByRegion.executeQuery();
        
            while(rs.next()) { 
                int start = rs.getInt(1);
                int stop = rs.getInt(2);
                char str = rs.getString(3).charAt(0);
                int poverlap = rs.getInt(4);
                
                ChipPetDatum intv = new ChipPetDatum(r.getGenome(), r.getChrom(), start, stop, str, expt, poverlap);
                intvs.addLast(intv);
            }
        }
        
        return intvs;
    }

    public void close() {
        try {
            loadByRegion.close(); loadByRegion = null;
        } catch (SQLException e) {
            e.printStackTrace();
        }
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }
    
}