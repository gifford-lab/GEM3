package edu.mit.csail.cgs.projects.readdb;

public class PairedHit {

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
        this.weight = weight;
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
                    Math.abs(weight - other.weight) < .001);
        } else {
            return false;
        }
    }   
    
    @Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + leftChrom;
		result = prime * result + leftPos;
		result = prime * result + (leftStrand ? 1231 : 1237);
		result = prime * result + rightChrom;
		result = prime * result + rightPos;
		result = prime * result + (rightStrand ? 1231 : 1237);
		return result;
	}

	public void flipSides() {
        int x = leftChrom;
        leftChrom = rightChrom;
        rightChrom = x;

        x = leftPos;
        leftPos = rightPos;
        rightPos = x;

        boolean b = leftStrand;
        leftStrand = rightStrand;
        rightStrand = b;

        short s = leftLength;
        leftLength = rightLength;
        rightLength = s;
    }
    public String toString() {
        return String.format("%d:%d,%d:%c and %d:%d,%d:%c weight %.2f",
                             leftChrom, leftPos, leftLength, leftStrand ? '+' : '-',
                             rightChrom, rightPos, rightLength, rightStrand ? '+' : '-',
                             weight);
    }

}