package edu.mit.csail.cgs.tools.chippet;

import java.io.PrintStream;
import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Interval;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLEntry;

public class ChromosomalBlatSummary {
    
    public static final int maxTreeSize = 1000;
	
	private Genome genome;
	private String chrom;
	private int chromLength;
    private String strand;
	private IntervalTree<String> intervalTree;
	private Map<String,Integer> keyCounts;
	
	public ChromosomalBlatSummary(Genome g, String c, String str) { 
		genome = g;
		chrom = c;
        strand = str;
		chromLength = genome.getChromLength(chrom);
		intervalTree = new IntervalTree<String>(new Interval(0, chromLength));
        keyCounts = new HashMap<String,Integer>();
	}

    public Genome getGenome() { return genome; }
    public String getChrom() { return chrom; }
    public String getStrand() { return strand; }
	public long getSize() { return intervalTree.getSize(); }
    public Iterator<Interval<String>> getIterator() { return intervalTree.getIterator(); }
    
    public void insertIntoDB(Connection cxn, int exptID) throws SQLException {  
		PreparedStatement cps = cxn.prepareStatement("select peakOverlap from chippetdata where " + 
				"expt=? and chromosome=? and startpos=? and stoppos=? and strand=?");
        PreparedStatement ps = cxn.prepareStatement("insert into chippetdata " +
                "(expt, chromosome, startpos, stoppos, strand, peakOverlap) values (?, ?, ?, ?, ?, ?)");
		PreparedStatement ups = cxn.prepareStatement("update chippetdata set peakOverlap=? where " + 
				"expt=? and chromosome=? and startpos=? and stoppos=? and strand=?");
        
        Iterator<Interval<String>> itr = intervalTree.getIterator();
        int chromID = genome.getChromID(chrom);
        
        long c = 0;
        
        while(itr.hasNext()) { 
            Interval<String> intv = itr.next();

			cps.setInt(1, exptID);
			cps.setInt(2, chromID);
			cps.setInt(3, intv.start);
			cps.setInt(4, intv.end);
            cps.setString(5, strand);

			ResultSet crs = cps.executeQuery();
			if(crs.next()) { 
				int po = crs.getInt(1);
				ups.setInt(1, po+1);
				ups.setInt(2, exptID);
				ups.setInt(3, chromID);
				ups.setInt(4, intv.start);
				ups.setInt(5, intv.end);
                ups.setString(6, strand);

				ups.executeUpdate();
			} else { 
				ps.setInt(1, exptID);
				ps.setInt(2, chromID);
				ps.setInt(3, intv.start);
				ps.setInt(4, intv.end);
                ps.setString(5, strand);
				ps.setInt(6, 1);
				
				ps.executeUpdate();
			}
			crs.close();

            c += 1;
        }
        
        ps.close();
		cps.close();
		ups.close();
        System.out.println("Inserted " + c + " values (" + intervalTree.getSize() + ")");
    }
    
    public void optimize() { 
        intervalTree.splitForDepth(10, maxTreeSize);
    }
    
    public void printSummary(PrintStream ps) { 
        ps.println("  Chrom Length: " + chromLength);
        intervalTree.printTreeSummary(1, ps);
    }
	
	public void addEntries(Collection<BlatPSLEntry> es) {
		for(BlatPSLEntry e : es) {
            addEntry(e);
		}
	}
	
	public void addEntry(BlatPSLEntry e) {
        Interval<String> intv = createInterval(e);
		intervalTree.addInterval(intv);
        String qname = intv.data;
        if(!keyCounts.containsKey(qname)) { keyCounts.put(qname, 0); }
        keyCounts.put(qname, keyCounts.get(qname)+1);
	}

    public int getCount(String key) { 
        return keyCounts.containsKey(key) ? keyCounts.get(key) : 0;
    }
    
    public int removeKeyedValue(String key) {
        int c = 0;
        if(keyCounts.containsKey(key)) {
            c = keyCounts.get(key);
            keyCounts.remove(key);
            intervalTree.removeValueWithData(key);
        }
        return c;
    }
    
    public boolean hasValue(String k) { return keyCounts.containsKey(k); }
    public Set<String> getKeys() { return keyCounts.keySet(); }
	
	public static Interval<String> createInterval(BlatPSLEntry es) {
		return new Interval<String>(es.getTstart(), es.getTend(), es.getQname());
	}

    public Collection<String> findOverlappingEntries(int start, int end) {
        LinkedList<Interval<String>> ints = new LinkedList<Interval<String>>();
        Interval query = new Interval(start, end);
        intervalTree.collectLeafIntervals(query, ints);
        
        LinkedList<String> entries = new LinkedList<String>();
        for(Interval<String> intv : ints) { 
            entries.addLast(intv.data);
        }
        
        return entries;
	}
}
