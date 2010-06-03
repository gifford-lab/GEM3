package edu.mit.csail.cgs.deepseq.features;

import java.util.List;

import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;

/**
 * A type of Peak found in CLIP-seq. Like every other Peak except stranded.
 * @author shaun
 *
 */
public class ClipSeqPeak extends EnrichedFeature{

	public ClipSeqPeak(){this('.');}
	public ClipSeqPeak(char s){this(null,s);}
	public ClipSeqPeak(Region c, char s){
		super(c);
		strand=s;
	}public ClipSeqPeak(Region c, double ip, double back, double z, double over, char str){
		super(c, ip, back, z, over, str);
	}
	
	public void setStrand(char s){strand=s;}
	
	public void updateStrand(List<ReadHit> overlapping){
		double forward=0, total=0;
		for(ReadHit r : overlapping){
			if(r.getStrand()=='+'){forward++;}
			total++;
		}
		if(total==0){strand='.';}
		else if(forward/total>0.5){strand='+';}
		else{strand='-';}
	}
	
	/** The peak string.  Fields are
     * - coordinates:strand
     * - region size
     * - location string
     * {- strand if CLIP}
     * - distance from location string to start of region
     * - number of hits (IP)
     * - number of hits (background)
     * - score (compares probability of read being in this region in IP and background)
     * - total hits (IP)
     * - total hits (background)
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
		int r2peak = peak.getLocation()-coords.getStart();
		return new String(coords.getLocationString()+":"+strand+"\t"+coords.getWidth()+"\t"+peak.getLocationString()+"\t"+r2peak+"\t"+String.format("%.1f", signalMaxHits)+"\t"+String.format("%.1f", backMaxHits)+"\t"+String.format("%.5e", score)+"\t"+String.format("%.1f", signalTotalHits)+"\t"+String.format("%.1f", backTotalHits)+"\t"+String.format("%.5f", overrep)+"\t"+gene+"\t"+distToGene+"\t" + annots.toString() + "\n");
	}
public String toGFF(){
		
		if (peak != null) {
			return new String(coords.getChrom()+"\tSEEDS\tpeak\t"+peak.getLocation()+"\t"+peak.getLocation()+"\t.\t"+strand+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalMaxHits+",Control="+backMaxHits);
		}else {
			return new String(coords.getChrom()+"\tSEEDS\tpeak\t"+peak.getLocation()+"\t"+peak.getLocation()+"\t.\t"+strand+"\t.\t"+"Note="+"Score:"+String.format("%.5e", score)+",Signal="+signalMaxHits+",Control="+backMaxHits);
		}
	}
	public String headString(){
		return new String("Region\tWidth\tPeak\tPeakOffset\tMaxSigHits\tMaxBackHits\tScore\tTotalSigHits\tTotalBackHits\tOverRep\tClosestGene\tTSSDist\tOtherAnnotations\n");
	}
	/**
	 * Override the Peak sequence maker to account for strandedness	 * 
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
		
		String seqName = new String(">"+peakWin.getLocationString()+":"+strand+"\t"+peakWin.getWidth());
		String sq = seqgen.execute(peakWin);
		if(strand=='-'){
			String tmp = SequenceUtils.reverseComplement(sq);
			sq=tmp;
		}
		seq = seqName +"\n"+ sq+"\n";
		return seq;
	}
}
