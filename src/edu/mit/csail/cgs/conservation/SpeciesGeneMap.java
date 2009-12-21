package edu.mit.csail.cgs.conservation;

/** 
 * @author tdanford
 *
 * An interface which returns GeneMap's for different species-pairs.
 */
public interface SpeciesGeneMap {
	public boolean hasMap(String startSpecies, String targetSpecies);
	public GeneMap getMap(String startSpecies, String targetSpecies);
}
