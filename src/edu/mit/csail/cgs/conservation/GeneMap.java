/**
 * 
 */
package edu.mit.csail.cgs.conservation;

import java.io.PrintStream;
import java.util.*;

/**
 * @author tdanford
 *
 * An interface for gene-to-gene maps between species.
 *
 */
public interface GeneMap {
	public Set<String> getDomainIDs();
	public Set<String> getRangeIDs();
	public Set<String> mapID(String id);
	public Set<String> mapIDs(Collection<String> ids);
	
	public void outputPairs(PrintStream ps);
}
