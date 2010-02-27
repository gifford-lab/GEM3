package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.util.ArrayList;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class MACSPeakRegion extends Region{
	int summit;
	int tags;
	double pvalue;
	double fold_enrichment;
	double fdr;
	public MACSPeakRegion(Genome g, String chr, int start, int end, 
			int summit, int tags, double pvalue, double fold_enrichment, double fdr){
		super(g, chr.replaceFirst("chr", ""), start, end);
		this.summit = summit;
		this.tags = tags;
		this.pvalue = pvalue;
		this.fold_enrichment = fold_enrichment;
		this.fdr = fdr;
	}
	public int getSummitRel() {
		return summit;
	}
	public int getSummitAbs() {
		return super.getStart()+summit;
	}
	public int getTags() {
		return tags;
	}
	//Convert to Region
	public Region getRegion(){
		return(new Region(this.getGenome(), this.getChrom(), this.getStart(), this.getEnd()));
	}//Convert to Peak (Point)
	public Point getPeak(){
		return(new Point(this.getGenome(), this.getChrom(), this.getSummitAbs()));
	}
	/*
	 * -10*log10(pvalue)
	 */
	public double getPvalue() {
		return pvalue;
	}
	public double getFold_enrichment() {
		return fold_enrichment;
	}
	/*
	 * %FDR
	 */
	public double getFdr() {
		return fdr;
	}
	public void setSummit(int summit) {
		this.summit = summit;
	}
	public void setTags(int tags) {
		this.tags = tags;
	}
	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}
	public void setFold_enrichment(double fold_enrichment) {
		this.fold_enrichment = fold_enrichment;
	}
	public void setFdr(double fdr) {
		this.fdr = fdr;
	}
	public static List<MACSPeakRegion> filterByPValue(List<MACSPeakRegion> peaks, double p_low, double p_high){
		List<MACSPeakRegion> peaks_filtered = new ArrayList<MACSPeakRegion>();
		for (MACSPeakRegion p:peaks){
			if (p.getPvalue()>p_high || p.getPvalue()<p_low){
				continue;
			}
			peaks_filtered.add(p);
		}
		return peaks_filtered;
	}
	public static List<MACSPeakRegion> filterByTags(List<MACSPeakRegion> peaks, int tag_low, int tag_high){
		List<MACSPeakRegion> peaks_filtered = new ArrayList<MACSPeakRegion>();
		for (MACSPeakRegion p:peaks){
			if (p.getTags()>tag_high || p.getTags()<tag_low){
				continue;
			}
			peaks_filtered.add(p);
		}
		return peaks_filtered;
	}
	public int compareByPValue(MACSPeakRegion macs) {	// descending p-value
		double diff = pvalue-macs.getPvalue();
		return diff==0?0:(diff<0)?1:-1;
	}
	
	public static List<MACSPeakRegion> filterByWidth(List<MACSPeakRegion> peaks, int width_low, int width_high){
		List<MACSPeakRegion> peaks_filtered = new ArrayList<MACSPeakRegion>();
		for (MACSPeakRegion p:peaks){
			if (p.getWidth()>width_high || p.getWidth()<width_low){
				continue;
			}
			peaks_filtered.add(p);
		}
		return peaks_filtered;
	}
	
	public String toMACS(){
		return getChrom()+"\t"+getStart()+"\t"+getEnd()+"\t"+(getEnd()-getStart())+"\t"+getSummitRel()
			+"\t"+getTags()+"\t"+getPvalue()+"\t"+getFold_enrichment()+"\t"+getFdr();
	}
}
