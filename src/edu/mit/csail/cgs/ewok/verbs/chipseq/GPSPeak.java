package edu.mit.csail.cgs.ewok.verbs.chipseq;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class GPSPeak extends Point{
	double strength;
	double controlStrength;
	double qvalue;
	double pvalue;
	double IPvsEMP;
	double IPvsCTR;
	private boolean jointEvent;		// 1 for joint event, 0 for unary, etc
	String nearestGene;
	int distance;
	Point EM_position;

	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double controlStrength, double qvalue, double pvalue, double IPvsEMP, 
			int joint, String nearestGene, int distance){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = controlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.IPvsEMP = IPvsEMP;
		this.IPvsCTR = 0;
		this.jointEvent = joint==1;
		this.nearestGene = nearestGene;
		this.distance = distance;
		this.EM_position = this;
	}

	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double controlStrength, double qvalue, double pvalue, double IPvsEMP, double IPvsCTR, 
			int joint, String nearestGene, int distance){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = controlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.IPvsEMP = IPvsEMP;
		this.IPvsCTR = IPvsCTR;
		this.jointEvent = joint==1;
		this.nearestGene = nearestGene;
		this.distance = distance;
		this.EM_position = this;
	}
	
	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double ctrlStrength, double qvalue, double pvalue, double IPvsEMP, double IPvsCTR){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = ctrlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.IPvsEMP = IPvsEMP;
		this.IPvsCTR = IPvsCTR;
		this.EM_position = this;
	}
	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double ctrlStrength, double qvalue, double pvalue, double IPvsEMP){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = ctrlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.IPvsEMP = IPvsEMP;
		this.IPvsCTR = 0;
		this.EM_position = this;
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
		return qvalue;
	}
	public double getPvalue() {
		return pvalue;
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
		this.qvalue = qvalue;
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
		return toString()+"\t"+strength+"\t"+controlStrength+"\t"+qvalue+"\t"+pvalue
		+"\t"+IPvsEMP+"\t"+jointEvent+"\t"+nearestGene+"\t"+distance;
	}
	
	public static String toGPS_Header(){
		return "GPS Event\tIP\tControl\tQ_-lg10\tP_-lg10\t"+
		"Shape\tJoint\tNearestGene\tDistance";
	}
	
	public String toGPS_short(){
		double fold = controlStrength==0 ? 9999:strength/controlStrength;
		String out = String.format("%s\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t%.3f", 
				toString(), strength, controlStrength, fold, qvalue, pvalue, IPvsEMP);
		return out;
	}
	
	public static String toGPS_short_Header(){
		return "Event Location\tIP\tControl\tIP/Ctrl\tQ_-lg10\tP_-lg10\tShape";
	}
	
	public String toGPS_motifShort(){
	  return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
	  +controlStrength;//+"\t"+qvalue+"\t"+shape+"\t"+shapeZ;
	}
}
