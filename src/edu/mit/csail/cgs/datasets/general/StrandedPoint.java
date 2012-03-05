package edu.mit.csail.cgs.datasets.general;

import java.util.regex.Matcher;

import edu.mit.csail.cgs.datasets.species.Genome;

public class StrandedPoint extends Point implements Stranded {
	
	private char strand;

	public StrandedPoint(Genome g, String c, int start, char str) {
		super(g, c, start);
		strand = str;
	}

	public char getStrand() {
		return strand;
	}

    public int hashCode() { 
        int code = super.hashCode();
        code += strand; code *= 37;
        return code;
    }

    public boolean equals(Object o) { 
        if(!(o instanceof StrandedPoint)) { return false; }
        StrandedPoint sp = (StrandedPoint)o;
        if(strand != sp.strand) { return false; }
        return super.equals(sp);
    }
    
    public String getLocationString() {
    	if (strand==' ') {
    		return super.getLocationString();
    	} else {
    		return super.getLocationString()+":"+strand;
    	}
      }
    
    /**
     * parses the input String into a StrandedPoint. Understands abbreviates in the
     * coordinates such as k and m. <br>
     * The method accepts the form: <blockquote>
     * 
     * <pre>
     * chromosome:start:strand
     * </pre>
     * 
     * </blockquote>
     */
    public static StrandedPoint fromString(Genome genome, String input) throws NumberFormatException {
      String pieces[] = input.split(":");
      char strand = ' ';
      if (pieces.length == 3 && pieces[2].length() == 1) {
        strand = pieces[2].charAt(0);
        input = pieces[0] + ":" + pieces[1];
      }

      String trimmed = input.trim();
      
      Matcher pointmatcher = POINT_PATTERN.matcher(trimmed);

      StrandedPoint output = null;
      
      if (pointmatcher.find()) {
        if (pointmatcher.groupCount() != 2) {
          return null;
        }
        String chromStr = pointmatcher.group(1);
        String locStr = pointmatcher.group(2);
        if (chromStr.startsWith("chr")) {
          chromStr = chromStr.substring(3, chromStr.length());
        }
        locStr = locStr.replaceAll(",", "");
        output = new StrandedPoint(genome, chromStr, Region.stringToNum(locStr), strand);
      }

      return output;
    }
}
