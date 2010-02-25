/*
 * Created on Jan 7, 2008
 */
package edu.mit.csail.cgs.utils.io.parsing.alignment;

import java.util.*;
import java.util.regex.*;

public class ElandHit extends AlignmentRecord {
    
    public static enum Code { U0, U1, U2, R0, R1, R2, QC, NM };
    
    private String name;
    private String querySequence;
    private Code code;
    private String targetSequence;
    private boolean strand;
    private int coordinate;
    private int[] mismatches;
    private String[] mismatchLocations;
    
    public ElandHit(int lineNum, String line) throws ElandParsingException { 
        String[] a = line.split("\\s+");
        
        if(a.length < 3) { 
            throw new ElandParsingException(lineNum, line);
        }
        
        name = removeFastaNamePattern(a[0]);
        querySequence = a[1];
        code = Code.valueOf(a[2]);
        
        mismatches = null;
        
        if(!a[2].startsWith("QC")) { 
            mismatches = new int[3];
            for(int i = 0; i < 3; i++) { 
                mismatches[i] = Integer.parseInt(a[i+3]);
            }

            targetSequence = null;
            strand = false;
            coordinate = -1;
            mismatchLocations = null;

            if(a[2].startsWith("U")) { 
                targetSequence = removeChromPattern(removeFastaPattern(a[6]));
                coordinate = Integer.parseInt(a[7]);
                strand = a[8].equals("F");
                mismatchLocations = new String[a.length-9+1];
                for(int i = 9; i < a.length; i++) { 
                    mismatchLocations[i-9] = a[i];
                }
            }
        }
    }
    
    public String getName() { return name; }
    public String getQuerySequence() { return querySequence; }
    public String getTargetSequence() { return targetSequence; }
    public int getCoordinate() { return coordinate; }
    public boolean getStrand() { return strand; }
    public Code getCode() { return code; }
    public int getMismatches(int m) { return mismatches[m]; }
    public String getMismatchLocation(int m) { return mismatchLocations[m]; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ElandHit)) { return false; }
        ElandHit h = (ElandHit)o;
        return h.name.equals(name);
    }
    
    public int hashCode() { 
        int code = 17;
        code += name.hashCode(); code *= 37;
        return code;
    }
    
    public String toString() { 
        return String.format("%s (%s)", name, code.toString());
    }
    
    public static Pattern fastaPattern = Pattern.compile("^(.*)\\.fa$");
    public static Pattern chromPattern = Pattern.compile("^chr(.*)$");
    public static Pattern fastaNamePattern = Pattern.compile("^>(.*)$");
    
    public static String removeChromPattern(String v) { 
        return removePattern(chromPattern, v);
    }
    
    public static String removeFastaPattern(String v) { 
        return removePattern(fastaPattern, v);
    }
    
    public static String removeFastaNamePattern(String v) { 
        return removePattern(fastaNamePattern, v);
    }
    
    public static String removePattern(Pattern p, String v) { 
        Matcher m = p.matcher(v);
        if(!m.matches()) { 
            return v;
        }
        return m.group(1);
    }
}
