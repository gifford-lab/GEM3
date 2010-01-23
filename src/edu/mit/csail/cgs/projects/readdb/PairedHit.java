package edu.mit.csail.cgs.projects.readdb;

public class PairedHit implements {

    public int leftChrom, rightChrom;
    public int leftPos, rightPos;
    public float weight;
    public boolean leftStrand, rightStrand;
    public short leftLength, rightLength;

    public PairedHit(int leftchrom, int leftpos, boolean leftstrand, short leftlen, 
                     int rightchrom, int rightpos, boolean rightstrand, short rightlen,
                     float weight) {
        leftChrom = leftchrom;
        rightChrom = rightchrom;
        leftPos = leftpos;
        rightPos = rightpos;
        weight = weight;
        leftStrand = leftstrand;
        rightStrand = rightstrand;
        leftLength = leftlen;
        rightLength = rightlen;            
    }

    public boolean equals(Object o) {
        if (o instanceof PairedHit) {
            PairedHit other = (PairedHit)o;
            return (leftChrom == other.leftChrom &&
                    rightChrom == other.rightChrom &&
                    leftPos == other.leftPos &&
                    rightPos == other.rightPos &&
                    leftStrand == other.leftStrand &&
                    rightStrand == other.rightStrand &&
                    weight == other.weight);
        }
    }   

}