package edu.mit.csail.cgs.tools.chippet;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.*;
import java.util.regex.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLEntry;
import edu.mit.csail.cgs.utils.io.parsing.alignment.BlatPSLParser;

public class PSLHitFinder {

	public static void main(String[] args) { 
		Pattern regpatt = Pattern.compile("([\\w\\d]+):(\\d+)-(\\d+)"); 
        File dir = new File(args[0]);
		int matches = 34;
		int gap = 10;

        try {
            Genome g = Organism.findGenome("mm6");
            Matcher regmatcher = regpatt.matcher(args[1]);
            if(!regmatcher.matches()) { throw new IllegalArgumentException(args[1]); }
            Region r = new Region(g, regmatcher.group(1), 
            			Integer.parseInt(regmatcher.group(2)), 
            			Integer.parseInt(regmatcher.group(3)));
            System.out.println("Matching " + r.getLocationString());
            
            ChipPetBlatPredicate pred = new ChipPetBlatPredicate(matches, gap);
            File[] pslfiles = getPSLFiles(dir);
            
            for(int i = 0; i < pslfiles.length; i++) { 
        		BlatPSLParser parser = new BlatPSLParser();
        		Iterator<BlatPSLEntry> entries = parser.parse(pslfiles[i], pred);
        		//System.out.println(String.format("--- %s ---", pslfiles[i].getName()));
                
        		while(entries.hasNext()) { 
        			BlatPSLEntry entry = entries.next();
                    String key = entry.getQname();
                    
                    String chrom = entry.getTname();
                    if(chrom.startsWith("chr")) { chrom=chrom.substring(3, chrom.length()); }
                    
                    int start = entry.getTstart(), end = entry.getTend();
                    Region tr = new Region(g, chrom, start, end);
                    
                    if(tr.overlaps(r)) { 
                    	System.out.println(key + "\t" + tr.getLocationString() + "\t" + pslfiles[i].getName());
					}
        		}                    
            }
            
        } catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static File[] getPSLFiles(File directory) throws IOException {
		FilenameFilter filter = new FilenameFilter() {
			public boolean accept(File f, String n) {
				return n.toUpperCase().endsWith(".PSL");
			} 
			
		};
		String[] names = directory.list(filter);
		File[] files = new File[names.length];
		for(int i = 0; i < names.length; i++) { files[i] = new File(directory, names[i]); }
		return files;
	}

}
