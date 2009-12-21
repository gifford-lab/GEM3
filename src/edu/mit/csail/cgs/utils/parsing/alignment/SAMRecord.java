package edu.mit.csail.cgs.utils.parsing.alignment;

public class SAMRecord {

    public String qname, tname, cigar, matetname, seq, qual;
    public int flags, pos, mapq, matepos, insertsize;

    /* this is all for the optional fields at the end.  nfields
       describes the number of valid entries in types, tags, and values
    */
    public int nfields;
    public char[] types;
    public String[] tags, values;

    /**
     * Assumes the quality string in this record is Phred formatted
     * (values between ! and I) and converts to the equivalent
     * quality string in Solexa format (< to h).
     */
    public void convertToSolexaQuals() {
        char[] newq = new char[qual.length()];
        for (int i = 0; i < qual.length(); i++) {
            int c = qual.charAt(i) - '!';
            double p = Math.pow(10, c / -10);
            c = '@' + (char) (-10 * Math.log10(p / (1-p)));
            if (c > 'h') {
                c = 'h';
            }
            newq[i] = (char)c;
        }
        qual = new String(newq);
    }
    /**
     * Assumes the quality string in this record is Solexa formatted
     * (values between < and h) and converts the equivalent
     * quality string in Solexa format (! to I).
     */
    public void convertToPhredQuals() {
        char[] newq = new char[qual.length()];
        for (int i = 0; i < qual.length(); i++) {
            int c = qual.charAt(i) - '@';
            double t = Math.pow(10, c / -10);
            double p =  t / (1 + t);
            c = '!' + (char) (-10 * Math.log10(p));
            newq[i] = (char)c;
        }
        qual = new String(newq);
    }
}