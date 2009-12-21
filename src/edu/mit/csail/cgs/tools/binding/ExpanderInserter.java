/*
 * Created on Nov 3, 2006
 */
package edu.mit.csail.cgs.tools.binding;

import java.util.*;
import java.util.regex.*;
import java.sql.*;
import java.io.*;

import edu.mit.csail.cgs.conservation.CustomMSPBindingGenerator;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.binding.*;

import edu.mit.csail.cgs.ewok.verbs.*;
/**
 * @author tdanford
 */
public class ExpanderInserter {
    
    public static void main(String[] args) {
    	ArgParser ap = new ArgParser(args);
    	if(!ap.hasKey("genome") || !ap.hasKey("expt") || !ap.hasKey("version")) { 
    		System.err.println("USAGE: ExpanderInserter --genome <genome name> " +
    				"--expt <experiment name> --version <version name> " +
    				"[--scanVersion <scan version name>] ");
    		return;
    	}
    	
        String genomeName = ap.getKeyValue("genome");
        String exptName = ap.getKeyValue("expt"), exptVersion = ap.getKeyValue("version");
        String scanVersion = ap.hasKey("scanVersion") ? ap.getKeyValue("scanVersion") : null;
        
        try {
            BindingScanLoader loader = new BindingScanLoader();            
            Genome genome = Organism.findGenome(genomeName);
            
            MSPLocator loc = new MSPLocator(genome, exptName, exptVersion);
            double p3Cutoff = 0.001;
            double pCutoff = 0.001;
            double nCutoff = 0.1;
            double twoCutoff = 0.005;
            
            Map<String,String> params = new HashMap<String,String>();
            params.put("p3", String.valueOf(p3Cutoff));
            params.put("p", String.valueOf(pCutoff));
            params.put("n", String.valueOf(nCutoff));
            params.put("two", String.valueOf(twoCutoff));
            
            Set<Region> regions = new HashSet<Region>();
            String version = scanVersion != null ? scanVersion : exptName + "," + exptVersion;

            CustomMSPBindingGenerator exp = 
                new CustomMSPBindingGenerator(loc.createObject(), p3Cutoff, pCutoff, nCutoff, twoCutoff);
            String type = exp.getClass().getSimpleName();
            
            int exptType = BindingScan.getLocatorType(loc);
            int[] exptIDs = BindingScan.getLocatorIDs(loc, loader.getConnection());
            BindingScan scan = new BindingScan(genome, type, version); 
            loader.insertScan(scan);            
            loader.insertNewExpts(scan, exptIDs, exptType);
            loader.insertNewParams(scan, params);
            loader.insertNewRegions(scan, regions);

            ExpanderInserter inserter = new ExpanderInserter(scan, loader, exp);
            
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            String line;

            Pattern p = Pattern.compile("([\\w\\d]+):([\\w\\d]+):(\\d+)-(\\d+)");
			Pattern allp = Pattern.compile("all ([\\w\\d]+)");
            
            while((line = br.readLine()) != null) { 
                Matcher m = p.matcher(line.trim());
				Matcher allm = allp.matcher(line.trim());
                if(m.matches()) { 
                    String gname = m.group(1);
                    Genome rg = Organism.findGenome(gname);
                    
                    String chrom = m.group(2);
                    int start = Integer.parseInt(m.group(3));
                    int end = Integer.parseInt(m.group(4));
                    
                    Region r = new Region(rg, chrom, start, end);
                    inserter.scanRegion(r);
                    
                    System.out.println(r.getGenome().getVersion() + ":" + r.getLocationString());
                } else if(allm.matches()) { 
					String gname = allm.group(1);
                    Genome rg = Organism.findGenome(gname);

                    for(String chrom : rg.getChromList()) { 
                    	int chromLength = rg.getChromLength(chrom);
                    	Region r = new Region(rg, chrom, 0, chromLength);
                    	
                    	inserter.scanRegion(r);
                        System.out.println(r.getGenome().getVersion() + ":" + r.getLocationString());
                    }
				}
            }
            
            loader.close();
            
        } catch (NotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private Expander<Region,BindingExtent> caller;
    private BindingScanLoader loader;
    
    private int dbid;
    private Map<Integer,BindingScan> scanMap;

    public ExpanderInserter(BindingScan bs, BindingScanLoader ld, Expander<Region,BindingExtent> c) 
        throws SQLException, UnknownRoleException {

        if(bs.getDBID() == -1) { 
            throw new IllegalArgumentException();
        }
        
        caller = c;
        loader = ld;
        
        dbid = bs.getDBID();
        
        scanMap = new HashMap<Integer,BindingScan>();
        scanMap.put(bs.getGenome().getDBID(), bs);
    }
    
    public void scanGenome(Genome g) throws SQLException {
        java.sql.Connection c = loader.getConnection();
        c.setAutoCommit(false);
        
        System.out.print(String.format("Scanning Genome (%s)", g.getVersion()));
        System.out.flush();

    	for(String chrom : g.getChromList()) { 
    		int length = g.getChromLength(chrom);
    		Region r = new Region(g, chrom, 0, length);
    		scanRegion(r);
            
            System.out.print(String.format(" %s", chrom));
            System.out.flush();
    	}
        System.out.println();

    	c.commit();
        c.setAutoCommit(true);
    }
    
    public void scanRegions(Collection<Region> regions) throws SQLException {
        java.sql.Connection c = loader.getConnection();
        c.setAutoCommit(false);
        
        for(Region r : regions) { scanRegion(r); }
        
        c.commit();
        c.setAutoCommit(true);
    }
    
    public void scanRegion(Region r) throws SQLException { 
        Genome rg = r.getGenome();
        BindingScan scan = getScan(rg);
        loader.insertNewRegion(scan, r);
        
        Iterator<BindingExtent> events = caller.execute(r);
        while(events.hasNext()) { 
            BindingExtent evt = events.next();
            loader.insertEvent(scan, evt);
        }        
    }
    
    public BindingScan getScan(Genome g) throws SQLException { 
    	if(scanMap.containsKey(g.getDBID())) { return scanMap.get(g.getDBID()); }
    	BindingScan scan = loader.loadScan(g, dbid);
    	scanMap.put(g.getDBID(), scan);
    	return scan;
    }
}
