package edu.mit.csail.cgs.deepseq.features;

import java.util.Collection;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;

/**
 * Peak is the parent class for the types of features found by statistical peak-finders.
 * Contains a peak point and various statistical descriptors. 
 * 
 * @author shaun
 *
 */
public class EnrichedNamedRegion extends Feature implements Comparable<EnrichedNamedRegion>{

	public String name;
	public double signalTotalHits;
	public double signalFPKM;
	public double backTotalHits;
	public double backFPKM;
	public double score;
	public double overrep;
	public char strand='.';
	
	/**
	 * For this EnrichedFeature: <br>
	 * region = null             <br>
	 * peak = null               <br>
	 * signalMaxHits = 0         <br>
	 * backMaxHits = 0           <br>
	 * signalTotalHits = 0       <br>
	 * backTotalHits = 0         <br>
	 * score = 0                 <br>
	 * overrep = 0               <br>
	 * strand = '.'             
	 */
    public EnrichedNamedRegion(){this(null);}
    
    
    /**
     * For this EnrichedFeature: <br>
	 * region = c                <br>
	 * signalMaxHits = 0         <br>
	 * backMaxHits = 0           <br>
	 * signalTotalHits = 0       <br>
	 * backTotalHits = 0         <br>
	 * score = 0                 <br>
	 * overrep = 0               <br>
	 * strand = '.'             
     * @param c region of this <tt>EnrichedFeature</tt>
     */
	public EnrichedNamedRegion(NamedRegion c){
		super(c);
		name=c.getName();
		signalFPKM = 0; backFPKM=0;
		signalTotalHits=0; backTotalHits=0;
		score=0; overrep=0;
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, double, double, double, double, char)
	 */
 	public EnrichedNamedRegion(NamedRegion c, double ip, double back, double s, double over){this(c, ip, back, s, over, '.');}
	
 	/**
     * For this EnrichedFeature: <br>
	 * region = c                <br>
	 * signalMaxHits = ip         <br>
	 * backMaxHits = back           <br>
	 * signalTotalHits = ip       <br>
	 * backTotalHits = back         <br>
	 * score = s                 <br>
	 * overrep = over               <br>
	 * strand = str              
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param ip number of IP (signal) hits
	 * @param back number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
 	 */
	public EnrichedNamedRegion(NamedRegion c, double ip, double back, double s, double over, char str){
		super(c);
		name=c.getName();
		signalFPKM=ip; backFPKM=back;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	public EnrichedNamedRegion(NamedRegion c, double ip, double back, double fpkmIP, double fpkmBack, double s, double over){this(c, ip, back, fpkmIP, fpkmBack, s, over, '.');}
	public EnrichedNamedRegion(NamedRegion c, double ip, double back, double fpkmIP, double fpkmBack, double s, double over, char str){
		super(c);
		name=c.getName();
		signalFPKM=fpkmIP; backFPKM=fpkmBack;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	
	//Rank according to increasing p-value (score)
	public int compareTo(EnrichedNamedRegion p) {
		if(score<p.score){return(-1);}
		else if(score>p.score){return(1);}
		else{return(0);}
	}
	
	 
	/** The EnrichedFeature string.  Fields are
     * - coordinates
     * - region size
     * - location string
     * - distance from location string to start of region
     * - max number of hits (IP)
     * - max number of hits (background)
     * - score (compares probability of read being in this region in IP and background)
     * - total number of hits (IP)
     * - total number of hits (background)
     * - overrepresentation factor
     * - closest gene
     * - distance to gene
     * - extra annotations if queried
     */
	public String toString() {
		StringBuilder annots = new StringBuilder();
        if (annotations != null) {
            for (Region r : annotations) {
                annots.append(r.toString() + ",");
            }
        }
        if (annots.length() > 1) {
            annots.setLength(annots.length() - 1);
        }
  
		String gene = nearestGene == null ? "NONE" : nearestGene.getName();
		
	    return new String(name+"\t"+coords.getLocationString()+"\t"+coords.getWidth()+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", signalTotalHits)+"\t"+String.format("%.1f", backTotalHits)+"\t"+String.format("%.5f", overrep)+"\t"+String.format("%.1f", signalFPKM)+"\t"+String.format("%.1f", backFPKM)+"\n");
		
	}
	//GFF3 (fill in later)
	public String toGFF(){
		return new String("Error:GFF writer not complete for EnrichedNamedRegion\n");
		//return new String(name+"\t"+coords.getLocationString()+"\t"+coords.getWidth()+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", signalTotalHits)+"\t"+String.format("%.1f", backTotalHits)+"\t"+String.format("%.5f", overrep)+"\t"+String.format("%.1f", signalFPKM)+"\t"+String.format("%.1f", backFPKM)+"\n");
	}
	
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of EnrichedFeature.toString().
	 */
	public String headString(){
		return new String("Name\tRegion\tWidth\tP-value\tSignalHits\tBackHits\tOverRep\tSignalFPKM\tBackFPKM\n");
	}
	
	
	/**
	 * Returns the sequence window (length=win) centered on the peak
	 */
	public String toSequence(int win){
		String seq;
		SequenceGenerator seqgen = new SequenceGenerator();
		Region reg=coords;
		String seqName = new String(">"+name+"\t"+reg.getLocationString()+"\t"+reg.getWidth());
		String sq = seqgen.execute(reg);
		seq = seqName +"\n"+ sq+"\n";
		return seq;
	}
}
