/*
 * Created on Jan 11, 2008
 */
package edu.mit.csail.cgs.datasets.chipseq;

import java.sql.SQLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;
import java.util.Iterator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.NotFoundException;

/**
 * @author tdanford
 *
 * ChipSeqLocator names a subset of the available ChipSeq experiments (rows from the 
 * chipseqexpts table), which satisfy the following conditions:
 *  
 * (1) the chipseqexpts.name matches the given exptName in the locator
 * (2) only the alignments with the given alignName are indicated
 * (3) if the locator's 'reps' field is non-empty, then only those replicates whose
 *     names appear in that set are indicated.
 */
public class ChipSeqLocator implements Comparable<ChipSeqLocator> {

    private String exptName, alignName;
    
    // This field has a particular semantics -- if it contains any entries, then 
    // the locator designates *only* those replicates (of the named experiment) 
    // that have the given alignment name too.
    // On the other hand, if "reps" is empty, then the locator designates *all* 
    // available replicates for which the given alignment name is valid.
    private Set<String> reps;
    

    public ChipSeqLocator(String ename) {
        exptName = ename;
        reps = new TreeSet<String>();
    }

    public ChipSeqLocator(String ename, String aname) {
        exptName = ename;
        alignName = aname;
        reps = new TreeSet<String>();
    }
    
    public ChipSeqLocator(String ename, String rname, String aname) {
        exptName = ename;
        alignName = aname;
        reps = new TreeSet<String>();
        reps.add(rname);
    }

    public ChipSeqLocator(String ename, Collection<String> rnames, String aname) {
        exptName = ename;
        alignName = aname;
        reps = new TreeSet<String>(rnames);
    }
    
    public boolean isSubset(ChipSeqLocator loc) { 
        if(!exptName.equals(loc.exptName)) { return false; }
        if(!alignName.equals(loc.alignName)) { return false; }
        for(String rep : reps) { 
            if(!loc.reps.contains(rep)) { return false; }
        }
        return true;
    }
    
    public String getExptName() { return exptName; }
    public String getAlignName() { return alignName; }
    
    public Collection<String> getReplicates() { return reps; }
    
    public String getReplicateString() {
        if(reps.isEmpty()) { return "all"; }
        StringBuilder sb = new StringBuilder();
        for(String rep : reps) { 
            if(sb.length() > 0) { sb.append(","); }
            sb.append(rep);
        }
        return sb.toString();
    }
    
    public int hashCode() { 
        int code = 17;
        code += exptName.hashCode(); code *= 37;
        code += alignName.hashCode(); code *= 37;
        
        for(String rep : reps) { code += rep.hashCode(); code *= 37; }
        return code;
    }
    
    public String toString() { 
    	return String.format("%s (rep: %s, align: %s)", exptName,
                getReplicateString(), alignName); 
    }

    public int compareTo(ChipSeqLocator other) {
        int c = exptName.compareTo(other.exptName);
        if (c == 0) {
            c = alignName.compareTo(other.alignName);
            if (c == 0) {
                Iterator<String> mine = reps.iterator();
                Iterator<String> others = other.reps.iterator();
                while (mine.hasNext() && others.hasNext()) {
                    c = mine.next().compareTo(others.next());
                    if (c != 0) {
                        break;
                    }
                }
                if (c == 0) {
                    if (mine.hasNext()) {
                        c = 1;
                    } else {
                        c = -1;
                    }
                }
            }
        }
        return c;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ChipSeqLocator)) { return false; }
        ChipSeqLocator loc =(ChipSeqLocator)o;
        if(!exptName.equals(loc.exptName)) { return false; }
        if(reps.size() != loc.reps.size()) { return false; }
        for(String rep : reps) { 
            if(!loc.reps.contains(rep)) { return false; }
        }
        if(!alignName.equals(loc.alignName)) { return false; }
        return true;
    }
    
    public Collection<ChipSeqAlignment> loadAlignments(ChipSeqLoader loader, Genome genome) 
    	throws SQLException, NotFoundException {
    	
        LinkedList<ChipSeqAlignment> alignments = new LinkedList<ChipSeqAlignment>();

        if(getReplicates().isEmpty()) { 
        	Collection<ChipSeqExpt> expts = loader.loadExperiments(getExptName());
        	for(ChipSeqExpt expt : expts) { 
        		ChipSeqAlignment alignment = loader.loadAlignment(expt, getAlignName(), genome);
        		if(alignment != null) { 
        			alignments.add(alignment);
        		}
        	}
        } else {
        	for(String repName : getReplicates()) { 
        		ChipSeqExpt expt = loader.loadExperiment(getExptName(), repName);
        		ChipSeqAlignment alignment = loader.loadAlignment(expt, getAlignName(), genome);
        		if(alignment != null) { 
        			alignments.add(alignment);
        		}
        	}
        }

        return alignments;
    }
}
