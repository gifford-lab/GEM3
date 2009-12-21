package edu.mit.csail.cgs.conservation;

import java.sql.SQLException;
import java.util.*;
import java.io.*;
import java.text.*;

import edu.mit.csail.cgs.datasets.function.DatabaseFunctionLoader;
import edu.mit.csail.cgs.datasets.function.FunctionLoader;
import edu.mit.csail.cgs.datasets.function.FunctionalUtils;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.probability.Hypergeometric;

public class BindingPartition implements edu.mit.csail.cgs.utils.Closeable {
    
    protected Hypergeometric hypgeom;
    protected SetTools<String> tools;
    protected NumberFormat nf;
    protected DecimalFormat sci;
	
    private String shortTitle, longTitle;
    private Vector<String> tags;
    private Map<String,Set<String>> blocks;
    private Set<String> total;
    
    private String fVersion;
    private FunctionLoader floader;
    private FunctionalUtils futils;
	
    public BindingPartition(String st, String lt, String fversion) {
        shortTitle = st; longTitle = lt;
        hypgeom = new Hypergeometric();
        tools = new SetTools<String>();
        blocks = new HashMap<String,Set<String>>();
        total = new HashSet<String>();
        tags = new Vector<String>();
        nf = new DecimalFormat("0.0000E0");
        sci = new DecimalFormat("0.0000E0");
        
        try {
            fVersion = fversion;
            floader = new DatabaseFunctionLoader();
            futils = new FunctionalUtils(floader, fversion);
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }
    public BindingPartition(String st, String lt, FunctionLoader loader, String fversion) {
        shortTitle = st; longTitle = lt;
        hypgeom = new Hypergeometric();
        tools = new SetTools<String>();
        blocks = new HashMap<String,Set<String>>();
        total = new HashSet<String>();
        tags = new Vector<String>();
        nf = new DecimalFormat("0.0000E0");
        sci = new DecimalFormat("0.0000E0");
        
        try {
            floader = loader;
            futils = new FunctionalUtils(floader, fversion);
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }

    public BindingPartition translate(Map<String,Set<String>> translation) { 
        BindingPartition bp = new BindingPartition(shortTitle, longTitle, fVersion);
        for(String tag : tags) { 
            HashSet<String> newBlock = new HashSet<String>();
            for(String id : blocks.get(tag)) { 
                if(translation.containsKey(id)) { 
                    newBlock.addAll(translation.get(id));
                }
            }
            bp.addBlock(tag, newBlock);
        }
        return bp;
    }
    
    public FunctionLoader getLoader() { return floader; }
    
    public void addBlock(String tag, Set<String> values) { 
        total.addAll(values);
        int removed = 0;
        for(String t : blocks.keySet()) {
            removed += tools.intersection(blocks.get(t), values).size();
            blocks.put(t, tools.subtract(blocks.get(t), values));
        }
        blocks.put(tag, new HashSet<String>(values));
        tags.add(tag);
        if(removed > 0) { 
            System.err.println("Removed " + removed + " values on tag \"" + tag + "\""); 
        }
    }
    
    public Enrichment getOverlapEnrichment(Set<String> sample, Set<String> bg) { 
    	Set<String> sampleids = new HashSet<String>();
    	Set<String> bgids = new HashSet<String>();
    	for(String tag : sample) { sampleids.addAll(getBlock(tag)); }
    	for(String tag : bg) { bgids.addAll(getBlock(tag)); }
    	
    	int N = total.size();
    	int n = sampleids.size();
    	int theta = bgids.size();
    	int x = tools.intersection(sampleids, bgids).size();
    	double log_pv = hypgeom.log_hypgeomPValue(N, theta, n, x);
    	
    	Enrichment e = new Enrichment(shortTitle, N, theta, n, x, log_pv);
    	return e;
    }
	
    public Set<String> getTotal() { return total; }
    public Vector<String> getBlockTags() { return tags; }
    public Set<String> getBlock(String tag) { return blocks.get(tag); }
    public String getShortTitle() { return shortTitle; }
    public String getLongTitle() { return longTitle; }
    public String getTitle() { return longTitle; }
    
    public String findTag(String id) { 
        for(String tag : tags) { 
            if(blocks.get(tag).contains(id)) { return tag; }
        }
        return null;
    }    
    
    public Enrichment[] getFDREnrichedGOCategories(String tag, double log_pvalue) throws SQLException { 
        return getFDREnrichedGOCategories(blocks.get(tag), log_pvalue);
    }
    
    public Enrichment[] getFDREnrichedGOCategories(Collection<String> tags, double log_pvalue) throws SQLException {
        HashSet<String> internal = new HashSet<String>();
        for(String t : tags) { 
            internal.addAll(blocks.get(t));
        }
        return getFDREnrichedGOCategories(internal, log_pvalue);
    }
    
    private Enrichment[] getFDREnrichedGOCategories(Set<String> internal, double log_fdr) throws SQLException { 
        Map<String,Enrichment> enrich = futils.calculateTotalEnrichments(total, internal);
        TreeSet<Enrichment> sorted = new TreeSet<Enrichment>(enrich.values());
        int count = 0;
        double factor = 0.0;
        int i = 0;
        int totalCount = sorted.size();
        
        if(totalCount > 0) {
            boolean keepTesting = true;
            for(Enrichment e : sorted) { 
                i += 1;
                if(keepTesting) { 
                    factor = (double)i / (double)totalCount;
                    double log_thresh = Math.log(factor) + log_fdr;
                    if(e.getLogPValue() <= log_thresh) { 
                        count += 1;
                    } else { 
                        keepTesting = false;
                    }
                }
            }
        }
       
        Enrichment[] sarray = sorted.toArray(new Enrichment[sorted.size()]);
        Enrichment[] array = new Enrichment[count];
        for(i = 0; i < count; i++) { 
            array[i] = sarray[i];
        }
        return array;
    }
    public void close() {
        futils.close();
        futils = null;
    }
    public boolean isClosed() {
        return (futils == null);
    }
    
    public Enrichment[] getEnrichedGOCategories(String tag, double log_pvalue) throws SQLException { 
        return getEnrichedGOCategories(blocks.get(tag), log_pvalue);
    }
    
    public Enrichment[] getEnrichedGOCategories(Collection<String> tags, double log_pvalue) throws SQLException {
        HashSet<String> internal = new HashSet<String>();
        for(String t : tags) { 
            internal.addAll(blocks.get(t));
        }
        return getEnrichedGOCategories(internal, log_pvalue);
    }
    
    private Enrichment[] getEnrichedGOCategories(Set<String> internal, double log_thresh) throws SQLException { 
        Map<String,Enrichment> enrich = futils.calculateTotalEnrichments(total, internal);
        TreeSet<Enrichment> sorted = new TreeSet<Enrichment>(enrich.values());
        TreeSet<Enrichment> threshed = new TreeSet<Enrichment>();
        
        int count = 0;
        double factor = 0.0;
        int i = 0;
        int totalCount = sorted.size();
        
        if(totalCount > 0) {
            boolean keepTesting = true;
            for(Enrichment e : sorted) { 
                i += 1;
                if(keepTesting) { 
                    factor = (double)i / (double)totalCount;
                    if(e.getLogPValue() <= log_thresh) { 
                        threshed.add(e);
                    }
                }
            }
        }
       
        Enrichment[] sarray = threshed.toArray(new Enrichment[threshed.size()]);
        return sarray;
    }
}
