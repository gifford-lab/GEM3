/*
 * Created on Mar 11, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.motifs;

import edu.mit.csail.cgs.ewok.verbs.Mapper;

/**
 * @author tdanford
 */
public class ReverseComplement implements Mapper<String,String> {
    
    private Mapper<Character,Character> complement;

    public ReverseComplement() {
        complement = new DNAComplement();
    }
    
    public ReverseComplement(Mapper<Character,Character> c) { 
        complement = c;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public String execute(String str) {
        StringBuilder sb = new StringBuilder();
        for(int i = str.length()-1; i >= 0; i--) { 
            sb.append(complement.execute(str.charAt(i)));
        }
        return sb.toString();
    }

    public static class DNAComplement implements Mapper<Character,Character> {
        
        public DNAComplement() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public Character execute(Character a) {
            switch(a.charValue()) { 
            case 'a':
                return 't';
            case 'A':
                return 'T';
            case 't':
                return 'a';
            case 'T':
                return 'A';
            case 'g':
                return 'c';
            case 'G':
                return 'C';
            case 'c':
                return 'g';
            case 'C':
                return 'G';
            }
            return 'N';
        }
        
    }

}
