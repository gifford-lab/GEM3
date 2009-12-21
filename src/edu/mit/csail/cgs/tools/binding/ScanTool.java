/*
 * Created on Nov 6, 2006
 */
package edu.mit.csail.cgs.tools.binding;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.conservation.CustomMSPBindingGenerator;
import edu.mit.csail.cgs.conservation.SimpleMSPBindingGenerator;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.utils.EwokBayes;
import edu.mit.csail.cgs.ewok.utils.EwokRosetta;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.CallerMapper;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author tdanford
 */
public class ScanTool {
	
    public static void main(String[] args) { 
        ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome") || !ap.hasKey("type")) { 
            System.err.println("Usage:\n" +
                               "ScanTool " +
                               "--species <organism name> " +
                               "--genome <genome version> " +
                               "--type [jbd|msp|hmm] " +
                               "--expt <db expt name> " +
                               "--version <db expt version> " + 
                               "--params [parameters]");
            return;
        }
        
        String species = ap.getKeyValue("species");
        String gname = ap.getKeyValue("genome");
        String callType = ap.getKeyValue("type");
        String params = ap.hasKey("params") ? ap.getKeyValue("params") : "";

        try {
            BindingScanLoader loader = new BindingScanLoader();
            Organism org = Organism.getOrganism(species);
            Genome genome = org.getGenome(gname);
            
            if(callType.equals("jbd")) { 
            	BindingParameters bps = new BindingParameters(params);
                CallerMapper cm = new ScanTool.JBDCallerMapper(bps);
                ScanTool scanner = new ScanTool(genome, loader, cm);
                if(params.length() > 0) { scanner.setParams(params); }
                
                String exptName = ap.getKeyValue("expt");
                String exptVersion = ap.getKeyValue("version");
                BayesLocator loc = new BayesLocator(genome, exptName, exptVersion);

                scanner.runScan(loc);
            } else if(callType.equals("msp")) { 
                CallerMapper cm = new ScanTool.MSPCallerMapper();
                ScanTool scanner = new ScanTool(genome, loader, cm);
                
                String exptName = ap.getKeyValue("expt");
                String exptVersion = ap.getKeyValue("version");
                MSPLocator loc = new MSPLocator(genome, exptName, exptVersion);

                scanner.runScan(loc);
            } else if(callType.equals("hmm")) { 
                BindingParameters bp = new BindingParameters(ap.getKeyValue("params"));
                CallerMapper cm = new ScanTool.HMMCallerMapper(bp);

                String exptName = ap.getKeyValue("expt");
                String exptVersion = ap.getKeyValue("version");
                ChipChipLocator loc = new ChipChipLocator(genome, exptName, exptVersion);

                ScanTool tool = new ScanTool(genome, loader, cm);
                tool.runScan(loc);
            }
			
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    private Genome genome;
    private BindingScanLoader loader;
    private java.sql.Connection cxn;
    private CallerMapper callerMapper;
    private String params;

    public ScanTool(Genome g,
                    BindingScanLoader l, 
                    CallerMapper cm) {
        genome = g;
        loader = l;
        cxn = loader.getConnection();
        callerMapper = cm;
        params = null;
    }
    
    public void setParams(String p) { 
    	params = p.replace(';', ',');
    }
    
    public void runScan(ExptLocator loc) throws SQLException, UnknownRoleException  { 
        LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
        locs.addLast(loc);
        runScan(locs);
    }

    public void runScan(Collection<? extends ExptLocator> locs) 
    	throws SQLException, UnknownRoleException {
    	
    	for(ExptLocator loc : locs) { 
            Expander<Region,BindingExtent> exp = callerMapper.execute(loc);
            String type = exp.getClass().getSimpleName();
            int exptType = BindingScan.getLocatorType(loc);
            int[] exptArray = BindingScan.getLocatorIDs(loc, cxn);
            Map<String,String> params = new HashMap<String,String>(callerMapper.getLastParams());
            String version = loc.getNameVersion().name + "," + loc.getNameVersion().version;
            if(params != null) { version += "," + params; }
    		
            BindingScan scan = new BindingScan(genome, type, version);
            loader.insertScan(scan);
            loader.insertNewExpts(scan, exptArray, exptType);
            loader.insertNewParams(scan, params);
    		
            ExpanderInserter expInsert = new ExpanderInserter(scan, loader, exp);
            expInsert.scanGenome(genome);
    		
            System.out.println("Scanned: " + loc.toString());
    	}
    }
    
    public static class HMMCallerMapper implements CallerMapper {
        
        private Map<String,String> params;
        private String config;
        
        public HMMCallerMapper(BindingParameters bp) { 
            params = new HashMap<String,String>();
            config = bp.getString("config", "edu.mit.csail.cgs.projects.ppg.3state_optimized_HMM");
            params.put("config", config);
        }

        public Map<String, String> getLastParams() {
            return params;
        }

        public Expander<Region, BindingExtent> execute(ExptLocator a) {
            ChipChipLocator loc = (ChipChipLocator)a;
            return new HMMDomainGenerator.BindingEventWrapper(new HMMDomainGenerator(config, loc));
        }
        
    }
    
    public static class DomainCallerMapper implements CallerMapper { 
        
        private Map<String,String> params;
        private CallerMapper internalCaller;
        private double density;
        private int minSize;
        
        public DomainCallerMapper(BindingParameters bp, CallerMapper cm) { 
            internalCaller = cm;
            params = new HashMap<String,String>();
            
            minSize = bp.getInt("size", 200);
            density = bp.getDouble("density", 1000.0);
            
            params.put("size", String.valueOf(minSize));
            params.put("density", String.valueOf(density));
        }
        
        public Map<String, String> getLastParams() {
            return params;
        }

        public Expander<Region, BindingExtent> execute(ExptLocator a) {
            Expander<Region,BindingExtent> gen = internalCaller.execute(a);
            BindingDomainGenerator domgen = new BindingDomainGenerator(gen, minSize, density);
            return domgen;
        }
    }
    
    public static class JBDCallerMapper implements CallerMapper {
    	
    	private Map<String,String> params;
    	private double probThresh, sizeThresh;
    	
    	public JBDCallerMapper() { 
            params = new HashMap<String,String>();
            probThresh = 0.2;
            sizeThresh = 2.0;
    		
            params.put("probThresh", String.valueOf(probThresh));
            params.put("sizeThresh", String.valueOf(sizeThresh));
    	}
    	
    	public JBDCallerMapper(BindingParameters bp) { 
            params = new HashMap<String,String>();
            probThresh = bp.getDouble("prob", 0.2);
            sizeThresh = bp.getDouble("str", 2.0);
    		
            params.put("probThresh", String.valueOf(probThresh));
            params.put("sizeThresh", String.valueOf(sizeThresh));    		
    	}

        public Map<String, String> getLastParams() {
            return new HashMap<String,String>(params);
        }

        public Expander<Region,BindingExtent> execute(ExptLocator a) {
            if(!(a instanceof BayesLocator)) { throw new IllegalArgumentException(); }
            BayesLocator loc = (BayesLocator)a;
            return new BayesBindingGenerator(loc.createObject(), probThresh, sizeThresh, true);
        } 
    	
    }
    
    public static class ChipChipCallerMapper implements CallerMapper {
    	
        private double outerThresh, peakThresh;
    	private Map<String,String> params;
    	
    	public ChipChipCallerMapper() { 
            outerThresh = 1.0;
            peakThresh = 2.0;
            params = new HashMap<String,String>();
            params.put("peakThresh", String.valueOf(peakThresh));
            params.put("outerThresh", String.valueOf(outerThresh));
    	}
    	
    	public ChipChipCallerMapper(BindingParameters bp) { 
            peakThresh = bp.getDouble("peakThresh", 2.0);
            outerThresh = bp.getDouble("outerThresh", 1.0);
            params = new HashMap<String,String>();
            params.put("peakThresh", String.valueOf(peakThresh));
            params.put("outerThresh", String.valueOf(outerThresh));    	    
    	}

		public Map<String, String> getLastParams() {
			return params;
		}

		public Expander<Region,BindingExtent> execute(ExptLocator a) {
            if(!(a instanceof ChipChipLocator)) { throw new IllegalArgumentException(); }
            ChipChipLocator al = (ChipChipLocator)a;
            return new ChipChipBindingGenerator(al.createObject(), outerThresh, peakThresh, true);
		} 
    	
    }
    
    public static class MSPCallerMapper implements CallerMapper {
        
        private double p3Cutoff, pCutoff, nCutoff, twoCutoff;
        private Map<String,String> params;
        
        public MSPCallerMapper() { 
            p3Cutoff = 0.001;
            pCutoff = 0.001;
            nCutoff = 0.1;
            twoCutoff = 0.005;            
            params = new HashMap<String,String>();
            params.put("p3", String.valueOf(p3Cutoff));
            params.put("p", String.valueOf(pCutoff));
            params.put("n", String.valueOf(nCutoff));
            params.put("two", String.valueOf(twoCutoff));
        }
        
        public MSPCallerMapper(BindingParameters bp) { 
            p3Cutoff = bp.getDouble("p3", 0.001);
            pCutoff = bp.getDouble("p", 0.001);
            nCutoff = bp.getDouble("n", 0.1);
            twoCutoff = bp.getDouble("two", 0.005);
            
            params = new HashMap<String,String>();
            params.put("p3", String.valueOf(p3Cutoff));
            params.put("p", String.valueOf(pCutoff));
            params.put("n", String.valueOf(nCutoff));
            params.put("two", String.valueOf(twoCutoff));        	
        }
        
        public MSPCallerMapper(double p3, double p, double n, double two) { 
            p3Cutoff = p3;
            pCutoff = p;
            nCutoff = n;
            twoCutoff = two;
            params = new HashMap<String,String>();
            params.put("p3", String.valueOf(p3Cutoff));
            params.put("p", String.valueOf(pCutoff));
            params.put("n", String.valueOf(nCutoff));
            params.put("two", String.valueOf(twoCutoff));
        }
        
        public Map<String,String> getLastParams() { 
            return params;
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Expander<Region, BindingExtent> execute(ExptLocator el) {
            if(!(el instanceof MSPLocator)) { throw new IllegalArgumentException(); }
            MSPLocator loc = (MSPLocator)el;
            
            CustomMSPBindingGenerator exp = 
                new CustomMSPBindingGenerator(loc.createObject(), p3Cutoff, pCutoff, nCutoff, twoCutoff);
            return exp;
        }            
        
    }

    public static class SimpleMSPCallerMapper implements CallerMapper {
        
        private double pCutoff;
        private Map<String,String> params;
        
        public SimpleMSPCallerMapper() { 
            pCutoff = 0.001;
            params = new HashMap<String,String>();
            params.put("p", String.valueOf(pCutoff));
        }
        
        public SimpleMSPCallerMapper(BindingParameters bp) { 
            pCutoff = bp.getDouble("p", 0.001);
            
            params = new HashMap<String,String>();
            params.put("p", String.valueOf(pCutoff));
        }
        
        public SimpleMSPCallerMapper(double p) { 
            pCutoff = p;
            params = new HashMap<String,String>();
            params.put("p", String.valueOf(pCutoff));
        }
        
        public Map<String,String> getLastParams() { 
            return params;
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Expander<Region, BindingExtent> execute(ExptLocator el) {
            if(!(el instanceof MSPLocator)) { throw new IllegalArgumentException(); }
            MSPLocator loc = (MSPLocator)el;
            
            SimpleMSPBindingGenerator exp = new SimpleMSPBindingGenerator(loc.createObject(), pCutoff);
            return exp;
        }            
        
    }
}
