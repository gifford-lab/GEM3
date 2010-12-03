package edu.mit.csail.cgs.deepseq.features;

import java.util.Collection;

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
public class EnrichedFeature extends Feature implements Comparable<EnrichedFeature>{
	public Point peak;
	public double signalMaxHits;
	public double backMaxHits;
	public double signalTotalHits;
	public double backTotalHits;
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
    public EnrichedFeature(){this(null);}
    
    
    /**
     * For this EnrichedFeature: <br>
	 * region = c                <br>
	 * peak = null               <br>
	 * signalMaxHits = 0         <br>
	 * backMaxHits = 0           <br>
	 * signalTotalHits = 0       <br>
	 * backTotalHits = 0         <br>
	 * score = 0                 <br>
	 * overrep = 0               <br>
	 * strand = '.'             
     * @param c region of this <tt>EnrichedFeature</tt>
     */
	public EnrichedFeature(Region c){
		super(c);
		peak=null;
		signalMaxHits=0;   backMaxHits=0;
		signalTotalHits=0; backTotalHits=0;
		score=0; overrep=0;
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, double, double, double, double, char)
	 */
 	public EnrichedFeature(Region c, double ip, double back, double s, double over){this(c, ip, back, s, over, '.');}
	
 	/**
     * For this EnrichedFeature: <br>
	 * region = c                <br>
	 * peak = start of region               <br>
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
	public EnrichedFeature(Region c, double ip, double back, double s, double over, char str){
		super(c);
		peak=new Point(c.getGenome(), c.getChrom(), c.getStart());
		signalMaxHits=ip;   backMaxHits=back;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * It uses the specified argument as a peak (instead of the start of the region).
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param peak_arg peak of this <tt>EnrichedFeature</tt>
	 * @param ip number of IP (signal) hits
	 * @param back number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
 	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, int peak_arg, double ip, double back, double s, double over, char str) {
		super(c);
		peak = new Point(c.getGenome(), c.getChrom(), peak_arg);
		signalMaxHits=ip; 	backMaxHits=back;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt> and
	 * makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, double ip, double back, double s, double over) {
		this(c, g, closestGeneDist, ip, back, s, over, '.');
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, double ip, double back, double s, double over, char str) {
		this(c, g, closestGeneDist, null, ip, back, s, over, str);
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>  
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, double ip, double back, double s, double over) {
		this(c, g, closestGeneDist, annots, ip, back, s, over, '.');
	}
	
	/**
	 * It stores information about the closest gene of this <tt>EnrichedFeature</tt>
	 * as well as the distance between them.          <br>
	 * Peak, here, is the start of the region.
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param g closest gene to this <tt>EnrichedFeature</tt>
	 * @param closestGeneDist distance between closest gene and this <tt>EnrichedFeature</tt>
	 * @param annots annotation of this <tt>EnrichedFeature</tt>
	 * @param ip number of IP (signal) hits
	 * @param back number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, double ip, double back, double s, double over, char str) {
		super(c, g, closestGeneDist, annots);
		peak=new Point(c.getGenome(), c.getChrom(), c.getStart());
		signalMaxHits=ip; 	backMaxHits=back;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt> and
	 * makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, int peak_arg, double ip, double back, double s, double over) {
		this(c, g, closestGeneDist, peak_arg, ip, back, s, over, '.');
	}

	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, int peak_arg, double ip, double back, double s, double over, char str) {
		this(c, g, closestGeneDist, null, peak_arg, ip, back, s, over, str);
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, int peak_arg, double ip, double back, double s, double over) {
		this(c, g, closestGeneDist, annots, peak_arg, ip, back, s, over, '.');
	}

	/**
	 * It stores information about the closest gene of this <tt>EnrichedFeature</tt>
	 * as well as the distance between them and uses the peak specified by the 
	 * corresponding argument.
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param g closest gene to this <tt>EnrichedFeature</tt>
	 * @param closestGeneDist distance between closest gene and this <tt>EnrichedFeature</tt>
	 * @param annots annotation of this <tt>EnrichedFeature</tt>
	 * @param peak_arg peak of this <tt>EnrichedFeature</tt>
	 * @param ip number of IP (signal) hits
	 * @param back number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, int peak_arg, double ip, double back, double s, double over, char str) {
		super(c, g, closestGeneDist, annots);
		peak = new Point(c.getGenome(), c.getChrom(), peak_arg);
		signalMaxHits=ip; 	backMaxHits=back;
		signalTotalHits=ip; backTotalHits=back;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * It uses the specified argument as a peak (instead of the start of the region).
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param peak_arg peak of this <tt>EnrichedFeature</tt>
	 * @param maxIP maximum number of IP (signal) hits
	 * @param maxBack maximum number of control (background) hits
	 * @param totalIP total number of IP (signal) hits
	 * @param totalBack total number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
 	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, int peak_arg, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over, char str) {
		super(c);
		peak = new Point(c.getGenome(), c.getChrom(), peak_arg);
		signalMaxHits=maxIP; 	backMaxHits=maxBack;
		signalTotalHits=totalIP; backTotalHits=totalBack;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt> and
	 * makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over) {
		this(c, g, closestGeneDist, maxIP, maxBack, totalIP, totalBack, s, over, '.');
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt>
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over, char str) {
		this(c, g, closestGeneDist, null, maxIP, maxBack, totalIP, totalBack, s, over, str);
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>  
	 * @see  edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over) {
		this(c, g, closestGeneDist, annots, maxIP, maxBack, totalIP, totalBack, s, over, '.');
	}
	
	/**
	 * It stores information about the closest gene of this <tt>EnrichedFeature</tt>
	 * as well as the distance between them.          <br>
	 * Peak, here, is the start of the region.
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param g closest gene to this <tt>EnrichedFeature</tt>
	 * @param closestGeneDist distance between closest gene and this <tt>EnrichedFeature</tt>
	 * @param annots annotation of this <tt>EnrichedFeature</tt>
	 * @param maxIP maximum number of IP (signal) hits
	 * @param maxBack maximum number of control (background) hits
	 * @param totalIP total number of IP (signal) hits
	 * @param totalBack total number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
	 */
	
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over, char str) {
		super(c, g, closestGeneDist, annots);
		peak=new Point(c.getGenome(), c.getChrom(), c.getStart());
		signalMaxHits=maxIP; 	backMaxHits=maxBack;
		signalTotalHits=totalIP; backTotalHits=totalBack;
		score=s; overrep=over;
		strand=str;
	}
	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt> and
	 * makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, int peak_arg, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over) {
		this(c, g, closestGeneDist, peak_arg, maxIP, maxBack, totalIP, totalBack, s, over, '.');
	}

	
	/**
	 * Sets the annotations of this <tt>EnrichedFeature</tt> to <tt>null</tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, int peak_arg, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over, char str) {
		this(c, g, closestGeneDist, null, peak_arg, maxIP, maxBack, totalIP, totalBack, s, over, str);
	}
	
	/**
	 * Makes use of the default value for the strand: <tt> strand = .  </tt>
	 * @see edu.mit.csail.cgs.deepseq.features.EnrichedFeature#EnrichedFeature(Region, Gene, int, Collection, int, double, double, double, double, double, double, char)
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, int peak_arg, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over) {
		this(c, g, closestGeneDist, annots, peak_arg, maxIP, maxBack, totalIP, totalBack, s, over, '.');
	}

	/**
	 * It stores information about the closest gene of this <tt>EnrichedFeature</tt>
	 * as well as the distance between them and uses the peak specified by the 
	 * corresponding argument.
	 * @param c region of this <tt>EnrichedFeature</tt>
	 * @param g closest gene to this <tt>EnrichedFeature</tt>
	 * @param closestGeneDist distance between closest gene and this <tt>EnrichedFeature</tt>
	 * @param annots annotation of this <tt>EnrichedFeature</tt>
	 * @param peak_arg peak of this <tt>EnrichedFeature</tt>
	 * @param maxIP maximum number of IP (signal) hits
	 * @param maxBack maximum number of control (background) hits
	 * @param totalIP total number of IP (signal) hits
	 * @param totalBack total number of control (background) hits
	 * @param s score (p-value) of this <tt>EnrichedFeature</tt>
	 * @param over over-representation of this <tt>EnrichedFeature</tt>
 	 * @param str strand of this <tt>EnrichedFeature</tt>
	 */
	public EnrichedFeature(Region c, Gene g, int closestGeneDist, Collection<Region> annots, int peak_arg, double maxIP, double maxBack, double totalIP, double totalBack, double s, double over, char str) {
		super(c, g, closestGeneDist, annots);
		peak = new Point(c.getGenome(), c.getChrom(), peak_arg);
		signalMaxHits=maxIP; 	backMaxHits=maxBack;
		signalTotalHits=totalIP; backTotalHits=totalBack;
		score=s; overrep=over;
		strand=str;
	}

	 
	//Rank according to increasing p-value (score)
	public int compareTo(EnrichedFeature p) {
		if(score<p.score){return(-1);}
		else if(score>p.score){return(1);}
		else{return(0);}
	}
	 public Point getPeak(){return(peak);}
	 
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
				
		if (peak != null) {
		  int r2peak = peak.getLocation()-coords.getStart();
	    return new String(coords.getLocationString()+"\t"+coords.getWidth()+"\t"+peak.getLocationString()+"\t"+r2peak+"\t"+String.format("%.1f", signalMaxHits)+"\t"+String.format("%.1f", backMaxHits)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", signalTotalHits)+"\t"+String.format("%.1f", backTotalHits)+"\t"+String.format("%.5f", overrep)+"\t"+gene+"\t"+distToGene+"\t" + annots.toString() + "\n");
		} else {
		  return new String(coords.getLocationString()+"\t"+coords.getWidth()+"\t"+String.format("%.1f", signalMaxHits)+"\t"+String.format("%.1f", backMaxHits)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", signalTotalHits)+"\t"+String.format("%.1f", backTotalHits)+"\t"+String.format("%.5f", overrep)+"\t"+gene+"\t"+distToGene+"\t" + annots.toString() + "\n");
		}
		
	}
	//GFF3 (fill in later)
	public String toGFF(){
		
		if (peak != null) {
			return new String(coords.getChrom()+"\tSEEDS\tpeak\t"+peak.getLocation()+"\t"+peak.getLocation()+"\t.\t.\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalMaxHits+",Control="+backMaxHits);
		}else {
			return new String(coords.getChrom()+"\tSEEDS\tpeak\t"+peak.getLocation()+"\t"+peak.getLocation()+"\t.\t.\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalMaxHits+",Control="+backMaxHits);
		}
	}
	
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of EnrichedFeature.toString().
	 */
	public String headString(){
		return new String("Region\tWidth\tPeak\tPeakOffset\tMaxSigHits\tMaxBackHits\tScore\tTotalSigHits\tTotalBackHits\tOverRep\tClosestGene\tTSSDist\tOtherAnnotations\n");
	}
	
	
	/**
	 * Returns the sequence window (length=win) centered on the peak
	 */
	public String toSequence(int win){
		String seq;
		SequenceGenerator seqgen = new SequenceGenerator();
		Region peakWin=null;
		if(win == -1){
			peakWin = coords;
		}else{
			int start = peak.getLocation()-((int)(win/2));
			if(start<1){start=1;}
			int end = peak.getLocation()+((int)win/2)-1;
			if(end>coords.getGenome().getChromLength(coords.getChrom())){end =coords.getGenome().getChromLength(coords.getChrom())-1;} 
			peakWin = new Region(coords.getGenome(), coords.getChrom(), start, end);
		}
		
		String seqName = new String(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth()+"\t"+peak.getLocation()+"\t"+new Integer(peak.getLocation()-coords.getStart()));
		String sq = seqgen.execute(peakWin);
		seq = seqName +"\n"+ sq+"\n";
		return seq;
	}
}
