/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.ewok.utils;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;

import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.probers.BayesImmediateProbeGenerator;
import edu.mit.csail.cgs.ewok.verbs.binding.*;
import edu.mit.csail.cgs.utils.*;

import edu.mit.csail.cgs.utils.database.*;

/**
 * @author tdanford
 */
public class EwokBayes extends EwokBase implements EwokBinding<BayesLocator,Probe> {
    
    public static void main(String[] args) { 
        ArgParser ap = new ArgParser(args);
        
        String orgName = "Homo sapiens";
        String genomeName = "hg18";
        if(ap.hasKey("organism")) { 
            orgName = ap.getKeyValue("organism");
        }
        if(ap.hasKey("genome")) { 
            genomeName = ap.getKeyValue("genome");
        }
		System.out.println("Organism: " + orgName);
		System.out.println("Genome: " + genomeName);

        EwokBayes et = new EwokBayes(orgName, genomeName);

        //binding_main(ap, et);
        try {
            Collection<BayesLocator> locs = et.getAllLocators();
            for(BayesLocator loc : locs) { 
                System.out.println("\"" + loc.name + "\" -- \"" + loc.version + "\"");
            }
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }
    
    public static void binding_main(ArgParser ap, EwokBayes et) { 
        String n = ap.getKeyValue("expt");
        String v = ap.getKeyValue("version");
        File input = new File(ap.getKeyValue("regions"));

        double pthresh = 0.5;
        double sthresh = 1.0;
        
        PeakCaller pc = et.getBayesPeakCaller(new BayesLocator(et.getGenome(), n, v), pthresh, sthresh);
        et.printRegionBindingAnnotation(input, pc);
    }
    
    public EwokBayes() {
        super();
    }
    
    public EwokBayes(String sp, String gn) {
        super(sp, gn);
    }
    
    public RegionProber<Probe> getProber(BayesLocator loc) { return getBayesProber(loc); }
    
    public Expander<Region,BindingExtent> getPeakCaller(BayesLocator loc, double[] params) {
        return new CustomBayesBindingGenerator(loc.createObject());
    }
    
    public RegionProber<Probe> getBayesProber(BayesLocator loc) { 
        BayesImmediateProbeGenerator gen = new BayesImmediateProbeGenerator(genome, loc);
        RegionProber<Probe> prober = new RegionProber.Wrapper<Probe>(gen);
        return prober;
    }
    
    public PeakCaller getBayesPeakCaller(BayesLocator loc, double pthresh, double sthresh) { 
        RegionProber prober = getBayesProber(loc);
        String n = loc.name, v = loc.version;
        BayesPeakFinder finder = new BayesPeakFinder(n, pthresh, sthresh);
        return new PeakCaller.FromFinder(prober, finder); 
    }
    
    public void printRegionBindingAnnotation(File regionFile, PeakCaller caller) {
        Iterator<Region> regions = super.getFileRegions(regionFile);
        
        BindingFilter bf = new BindingFilter(caller);
        Iterator<Pair<Region,Integer>> annRegions = 
            new MapperIterator<Region,Pair<Region,Integer>>(new FilterValueMapper<Region,Region>(bf), regions);
        
        while(annRegions.hasNext()) {
            Pair<Region,Integer> p = annRegions.next();
            System.out.println(p.getLast() + "\t" + p.getFirst());
        }
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.analysis.ewok.EwokBinding#getAllLocators()
     */
    public Collection<BayesLocator> getAllLocators() throws SQLException, UnknownRoleException {
        //int species = getGenome().getSpeciesDBID();
        int genomeID = getGenome().getDBID();
        
        LinkedList<BayesLocator> locs = new LinkedList<BayesLocator>();
        java.sql.Connection c = DatabaseFactory.getConnection(ExptLocator.dbRole);
        
        Statement s = c.createStatement();
        ResultSet rs = s.executeQuery("select ma.name, ma.version from bayesanalysis ma, bayesToGenome b2g " + 
                " where ma.id=b2g.analysis and b2g.genome=" + genomeID);

        while(rs.next()) { 
            String name2 = rs.getString(1);
            String version2 = rs.getString(2);
            BayesLocator locb = new BayesLocator(getGenome(), name2, version2);
            locs.addLast(locb);
        }

        rs.close();
        s.close();
        DatabaseFactory.freeConnection(c);
        
        return locs;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.analysis.ewok.EwokBinding#getBase()
     */
    public EwokBase getBase() {
        return this;
    }    
}
