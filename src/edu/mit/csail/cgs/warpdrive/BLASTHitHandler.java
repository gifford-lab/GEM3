/*
 * Created on Mar 1, 2007
 */
package edu.mit.csail.cgs.warpdrive;

import java.util.*;

import java.io.PrintStream;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.xml.XMLTreeElement;
import edu.mit.csail.cgs.utils.xml.XMLTreeHandler;
import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class BLASTHitHandler extends XMLTreeHandler {
    
    public static Set<String> tags;
    
    static { 
        tags = new HashSet<String>();
        tags.add("Hit");
        tags.add("Hsp");
        tags.add("Iteration");
        tags.add("Parameters");
        tags.add("Statistics");
    }
    
    public BLASTHitHandler() {
        super();
        for(String t : tags) { addTag(t); }
    }
    
    public LinkedList<ScoredStrandedRegion> extractHitRegions(Genome g) { 
        LinkedList<ScoredStrandedRegion> regions = new LinkedList<ScoredStrandedRegion>();
        
        LinkedList<XMLTreeElement> hits = getTree().collectElements("Hit");
        for(XMLTreeElement hit : hits) { 
            
            String chrom = hit.values.get("Hit_def");
            
            LinkedList<XMLTreeElement> hsps = hit.collectElements("Hsp");
            for(XMLTreeElement hsp : hsps) { 

                int from = Integer.parseInt(hsp.values.get("Hsp_hit-from"));
                int to = Integer.parseInt(hsp.values.get("Hsp_hit-to"));
                double score = Double.parseDouble(hsp.values.get("Hsp_score"));

                char strand = '+';
                
                if(from > to) { 
                    int temp = from;
                    from = to;
                    to = from;
                    strand = '-';
                }
                
                ScoredStrandedRegion r = new ScoredStrandedRegion(g, chrom, from, to, score, strand);
                regions.addLast(r);
            }
        }
        
        return regions;
    }    
}
