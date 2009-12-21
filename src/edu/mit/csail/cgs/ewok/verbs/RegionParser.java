/*
 * Created on Mar 9, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.util.regex.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class RegionParser implements Mapper<String,Region> {
    
    private static Pattern regPatt;
    
    static { 
        regPatt = Pattern.compile("(\\w+):(\\d+)-(\\d+)");
    }
    
    private Genome genome;
    private int chromIndex, startIndex, endIndex, nameIndex, minLength;

    public RegionParser(Genome g) {
        genome = g;
        chromIndex = 0;
        startIndex = 1;
        endIndex = 2;
        nameIndex = 3;
        minLength = (Math.max(chromIndex, Math.max(startIndex, endIndex))) + 1;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public Region execute(String input) {
        String[] array = input.split("\\s+");
        String chrom = array[chromIndex];
        
        Matcher m = regPatt.matcher(chrom);
        if(m.matches()) { 
            chrom = m.group(1);
            int start = Integer.parseInt(m.group(2));
            int end = Integer.parseInt(m.group(3));
            return new Region(genome, chrom, start, end);
        } else { 
            if(array.length >= minLength) { 
                int start = Integer.parseInt(array[startIndex]);
                int end = Integer.parseInt(array[endIndex]);
                if(nameIndex < array.length) {
                    return new NamedRegion(genome, chrom, start, end, array[nameIndex]);
                } else { 
                    return new Region(genome, chrom, start, end);
                }
            } else { 
                System.err.println("Line \"" + input + "\" doesn't have the correct length (" + minLength + ")");
                return null;
            }
        }
    }

}
