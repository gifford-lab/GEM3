package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class StrandedRegion extends Region implements Stranded {
    
    protected char strand;
    
    public StrandedRegion(StrandedRegion nr) { 
        super(nr);
        strand = nr.strand;
    }
    
    public StrandedRegion(Region r) { 
        super(r);
        strand = 1;
    }
    
    public StrandedRegion(Region r, char st) { 
    	super(r);
    	strand = st;
    }
    
    public StrandedRegion(Genome g, String c, int start, int end, char str) {
        super(g,c,start,end);
        strand = str;
    }

    public char getStrand() { return strand; }
    public int getFivePrime() {
        if (strand == '+') {
            return getStart();
        } else if (strand == '-') {
            return getEnd();
        } else {
            throw new IllegalArgumentException("Unknown strand " + getStrand());
        }
    }
    public int getThreePrime() {
        if (strand == '-') {
            return getStart();
        } else if (strand == '+') {
            return getEnd();
        } else {
            throw new IllegalArgumentException("Unknown strand " + getStrand());
        }
    }
    
    public String toString() { 
        return getLocationString() + ":" + strand;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof StrandedRegion)) { return false; }
        StrandedRegion nr = (StrandedRegion)o;
        if(nr.strand != strand) { return false; }
        return super.equals(nr);
    }

    public StrandedRegion expand(int upstream, int downstream) {
        if (strand == '+') {
            int ns = getStart() - upstream;
            int ne = getEnd() + downstream;
            if (ns < 1) {ns = 1;}
            return new StrandedRegion(getGenome(),getChrom(),ns,ne,strand);
        } else if (strand == '-') {
            int ns = getStart() - downstream;
            int ne = getEnd() + upstream;
            if (ns < 1) {ns = 1;}
            return new StrandedRegion(getGenome(),getChrom(),ns,ne,strand);                
        } else {
            throw new IllegalArgumentException("Strand isn't + or - so I don't know what to do");
        }

    }

    public static StrandedRegion fromString(Genome genome, String input) throws NumberFormatException {
        String trimmed = input.trim();
        String regionregex = "^\\s*([\\w\\d]+):\\s*([,\\d]+[mMkK]?)\\s*-\\s*([,\\d]+[mMkK]?):([\\+\\-])\\s*";
        Pattern regionpattern = Pattern.compile(regionregex);
        Matcher regionmatcher = regionpattern.matcher(trimmed);
        if (regionmatcher.find()) {
            if(regionmatcher.groupCount() != 4) { return null; }
            String chromStr = regionmatcher.group(1);
            String startStr = regionmatcher.group(2);
            String endStr = regionmatcher.group(3);
            char strand = regionmatcher.group(4).charAt(0);
            if(chromStr.startsWith("chr")) { chromStr = chromStr.substring(3, chromStr.length()); }
            startStr = startStr.replaceAll(",", "");
            endStr = endStr.replaceAll(",", "");
            int start = Math.min(stringToNum(startStr),
                                 stringToNum(endStr));
            int end = Math.max(stringToNum(startStr),
                                 stringToNum(endStr));
            return new StrandedRegion(genome,
                                      chromStr,
                                      start,
                                      end,
                                      strand);
        } 
        return null;
    }
    
    public int hashCode() { 
        int code = super.hashCode();
        code += strand; code *= 37;
        return code; 
    }
}
