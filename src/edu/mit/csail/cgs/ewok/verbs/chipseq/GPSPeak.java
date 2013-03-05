package edu.mit.csail.cgs.ewok.verbs.chipseq;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class GPSPeak extends Point{
	public static char UNKNOWN_STRAND = '*';
	double strength;
	double controlStrength;
	double expectedStrength;
	double qv_lg10;
	double pv_lg10;
	double pvalue;
	double pois_pv_lg10;
	double IPvsEMP;
	double IPvsCTR;
	private boolean jointEvent;		// 1 for joint event, 0 for unary, etc
	String nearestGene;
	int distance;
	Point EM_position;
	String kmer;
	int kmerGroupCount;
	char kmerStrand = UNKNOWN_STRAND;
	double kmerStrength;
	String boundSequence;
	public String getKmer() {return kmer;}
	public void setKmer(String km) {kmer=km;}
	public int getKmerGroupCount() {	return kmerGroupCount;}
	public double getKmerStrength() {return kmerStrength;}
	public char getKmerStrand() {return kmerStrand;}
	public String getBoundSequence() {return boundSequence;}
	public void setBoundSequence(String bs) {boundSequence=bs;}


	// GPS output format 2011-07-25	
	// Position	     IP	Control	   Fold	Expectd	Q_-lg10	P_-lg10	P_poiss	IPvsEMP	IPvsCTR	Kmer	Count	Strength	BoundSequence	EnrichedHGP
	public GPSPeak(Genome g, String chr, int pos, double ipStrength, double ctrlStrength, double expectedStrength,
			double qv_lg10, double pvalue, double pv_lg10, double pois_pvalue, double IPvsEMP, double IPvsCTR, 
			String kmer, int kmerCount, char kmerStrand, String boundSequence){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = ctrlStrength;
		this.expectedStrength = expectedStrength;
		this.qv_lg10 = qv_lg10;
		this.pvalue = pvalue;
		this.pv_lg10 = pv_lg10;
		this.pois_pv_lg10 = pois_pvalue;
		this.IPvsEMP = IPvsEMP;
		this.IPvsCTR = IPvsCTR;
		this.EM_position = this;
		this.kmer = kmer;
		this.kmerGroupCount = kmerCount;
		this.kmerStrand = kmerStrand;
		this.boundSequence = boundSequence;
	}
	
	public boolean isJointEvent() {
		return jointEvent;
	}

	public double getStrength() {
		return strength;
	}
	public double getControlStrength() {
		return controlStrength;
	}
	public Point getEM_position() {
		return EM_position;
	}
	public double getQvalue() {
		return qv_lg10;
	}
	public double getPvalue() {
		return pv_lg10;
	}
	public double getPV_lg10() {
		return pv_lg10;
	}
	public double getShape() {
		return IPvsEMP;
	}
	public double getIPvsCTR() {
		return IPvsCTR;
	}
	public void setStrength(double strength) {
		this.strength = strength;
	}
	public void setControlStrength(double controlStrength) {
		this.controlStrength = controlStrength;
	}
	public void setEM_position(Point em_position) {
		EM_position = em_position;
	}
	public void setQvalue(double qvalue) {
		this.qv_lg10 = qvalue;
	}
	public void setShape(double shape) {
		this.IPvsEMP = shape;
	}

	public String getNearestGene() {
		return nearestGene;
	}
	public int getDistance() {
		return distance;
	}
	public void setNearestGene(String nearestGene) {
		this.nearestGene = nearestGene;
	}
	public void setDistance(int distance) {
		this.distance = distance;
	}
	public int compareByPValue(GPSPeak p) {
		double diff = getPvalue()- p.getPvalue();
		return	diff==0?0:(diff<0)?-1:1;	//p-value: ascending
	}
	public int compareByIPStrength(GPSPeak p) {
		double diff = getStrength()- p.getStrength();
		return	diff==0?0:(diff<0)?1:-1;	//ip strength: descending
	}
	public String toGPS(){
		return toString()+"\t"+strength+"\t"+controlStrength+"\t"+qv_lg10+"\t"+pv_lg10
		+"\t"+IPvsEMP+"\t"+jointEvent+"\t"+nearestGene+"\t"+distance;
	}
	
	public static String toGPS_Header(){
		return "GPS Event\tIP\tControl\tQ_-lg10\tP_-lg10\t"+
		"Shape\tJoint\tNearestGene\tDistance";
	}
	
	public String toGPS_short(){
		double fold = controlStrength==0 ? 9999:strength/controlStrength;
		String out = String.format("%s\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t%.3f", 
				toString(), strength, controlStrength, fold, qv_lg10, pv_lg10, IPvsEMP);
		return out;
	}
	
	public static String toGPS_short_Header(){
		return "Event Location\tIP\tControl\tIP/Ctrl\tQ_-lg10\tP_-lg10\tShape";
	}
	
	public String toGPS_motifShort(){
	  return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
	  +controlStrength;//+"\t"+qvalue+"\t"+shape+"\t"+shapeZ;
	}
	/** 
	 * Output ENCODE narrow peak format<br>
	 * The score is based on q-value (corrected from the larger p-value from binomial p-value (IP vs Control) and Poisson p-value (IP only)
	 * @return
	 */
	public String toNarrowPeak(double score){
		return "chr"+getChrom()+"\t"+(getLocation()-100)+"\t"+(getLocation()+101)+"\tpeak\t"+score+"\t"+"."+"\t"+strength+"\t"+(pv_lg10>pois_pv_lg10?pois_pv_lg10:pv_lg10)+"\t"+qv_lg10
		+"\t"+100;
	}
}
