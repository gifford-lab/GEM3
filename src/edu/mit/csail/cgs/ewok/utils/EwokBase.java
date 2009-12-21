/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.ewok.utils;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.locators.BayesLocator;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.iterators.SingleIterator;

/**
 * @author tdanford
 */
public class EwokBase {
    
    protected Organism org;
    protected Genome genome;
    private RefGeneByNameGenerator nameExpander;
    
    public EwokBase() {
        try {
            org = Organism.getOrganism("Saccharomyces cerevisiae");
            genome = org.getGenome("sacCer1");
            nameExpander = new RefGeneByNameGenerator(genome);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public EwokBase(String sp, String gn) {
        try {
            org = Organism.getOrganism(sp);
            genome = org.getGenome(gn);
            nameExpander = new RefGeneByNameGenerator(genome);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public EwokBase(Genome g) { 
        try {
            org = Organism.getOrganism(g.getSpecies());
            genome = g;
            nameExpander = new RefGeneByNameGenerator(genome);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }        
    }
    
    public int hashCode() { 
        int code = 17;
        code += genome.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof EwokBase)) { return false; }
        EwokBase eb = (EwokBase)o;
        if(!genome.equals(eb.genome)) { return false; }
        return true;
    }
    
    public String toString() { return "EwokBase: " + org.getName() + " (" + genome.getVersion() + ")"; }
    
    public Organism getOrganism() { return org; }
    public Genome getGenome() { return genome; }
    
    public Iterator<NamedRegion> getChroms() { 
        Iterator<NamedRegion> chromItr = new ChromRegionIterator(genome);
        return chromItr;
    }    
    
    public Iterator<Region> getFileRegions(File regionFile) {
        Iterator<File> file = new SingleIterator<File>(regionFile);
        Iterator<String> lines = new ExpanderIterator<File,String>(new FileLineExpander(), file);
        Iterator<Region> regions = new MapperIterator<String,Region>(new RegionParser(genome), lines);
        return regions;
    }
    
    public Iterator<Gene> loadGenesByName(String n) {
        /*
        Iterator<String> names = new SingleIterator<String>(n);
        Iterator<Gene> genes = new ExpanderIterator<String,Gene>(new NameToGene(genome), names);
        return genes;
        */
        return loadGenesByName(new SingleIterator<String>(n));
    }
    
    public Iterator<Gene> loadGenesByName(Iterator<String> names) { 
        Iterator<Gene> genes = new ExpanderIterator<String,Gene>(nameExpander, names);
        return genes;        
    }
    
    public Iterator<Gene> loadGenesByID(String id) { 
        Iterator<Gene> baseGenes = loadGenesByName(id);
        Iterator<Gene> genes = new FilterIterator<Gene,Gene>(new GeneIDFilter(id), baseGenes);
        return genes;
    }
    
    public static class GeneSourceFilter implements Filter<Gene,Gene> { 
        private String source;
        
        public GeneSourceFilter(String s) { source = s; }
        
        public Gene execute(Gene a) { 
            if(a.getSource().equals(source)) { 
                return a;
            }
            return null;
        }
    }
    
    public static class GeneIDFilter implements Filter<Gene,Gene> { 
        private String id;
        public GeneIDFilter(String i) { id = i; }
        public Gene execute(Gene a) { 
            if(a.getID().equals(id)) { return a; }
            return null;
        }
    }
    
    public Iterator<Region> loadRegions(File f) { 
    	Iterator<File> file = new SingleIterator<File>(f);
    	Iterator<String> lines = new ExpanderIterator<File,String>(new FileLineExpander(), file);
    	RegionParser rp = new RegionParser(genome);
    	Iterator<Region> regions = new MapperIterator<String,Region>(rp, lines);
    	return regions;
    }
    
}
