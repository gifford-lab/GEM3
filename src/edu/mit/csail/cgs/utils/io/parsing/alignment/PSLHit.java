package edu.mit.csail.cgs.utils.io.parsing.alignment;

public class PSLHit {
	public int match, mismatch, repmatch, n, qgapcount, qgapbases, tgapcount, tgapbases;
	public char strand;
	public String qname, tname;
	public int qsize, qstart, qend, tsize, tstart, tend;
	public int[] blockSizes, qStarts, tStarts;

	public int matchSize() { return match + mismatch; }
	
	public double getMatchedIdentity() { 
		return (double)(match) / (double)(match+mismatch);
	}
	
	public double getQueryIdentity() { 
		return (double)(match) / (double)(qsize);
	}
    public String toString() {
        StringBuilder bsizes = new StringBuilder();
        StringBuilder qstarts = new StringBuilder();
        StringBuilder tstarts = new StringBuilder();
        for (int i = 0; i < blockSizes.length; i++) {
            bsizes.append(blockSizes[i]);
            qstarts.append(qStarts[i]);
            tstarts.append(tStarts[i]);
            if (i != blockSizes.length - 1) {
                bsizes.append(",");
                qstarts.append(",");
                tstarts.append(",");
            }
        }

        return String.format("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t" +
                             "%s\t%d\t%d\t%d" +
                             "%s\t%d\t%d\t%d" +
                             "%d\t%s\t%s\t%s",
                             match,
                             mismatch,
                             repmatch,
                             n,
                             qgapcount,
                             qgapbases,
                             tgapcount,
                             tgapbases,
                             strand,
                             qname,
                             qsize,
                             qstart,
                             qend,
                             tname,
                             tsize,
                             tstart,
                             tend,
                             blockSizes.length,
                             bsizes,
                             qstarts,
                             tstarts);
    }
}
