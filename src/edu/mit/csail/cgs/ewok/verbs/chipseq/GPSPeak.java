package edu.mit.csail.cgs.ewok.verbs.chipseq;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class GPSPeak extends Point{
	double strength;
	double controlStrength;
	Point EM_position;
	double qvalue;
	double pvalue;
	double shape;
	double shapeZ;
	String nearestGene;
	int distance;
	double mixProb;
	int unaryEvent;		// 1 for unary event, 0 for binary, etc
	
	 public GPSPeak(Genome g, String chr, int pos, int EM_pos, double strength, 
      double controlStrength, double qvalue, double shape, double shapeZ,
      double mixProb, String nearestGene, int distance){
    super(g, chr.replaceFirst("chr", ""), pos);
    EM_position = new Point(g, chr.replaceFirst("chr", ""), EM_pos);
    this.strength = strength;
    this.controlStrength = controlStrength;
    this.qvalue = qvalue;
    this.shape = shape;
    this.shapeZ = shapeZ;
    this.mixProb = mixProb;
    this.nearestGene = nearestGene;
    this.distance = distance;
    
    this.pvalue = Double.NaN;
  }
	
	public GPSPeak(Genome g, String chr, int pos, int EM_pos, double strength, 
			double controlStrength, double qvalue, double pvalue, double shape, double shapeZ,
			double mixProb, String nearestGene, int distance){
		super(g, chr.replaceFirst("chr", ""), pos);
		EM_position = new Point(g, chr.replaceFirst("chr", ""), EM_pos);
		this.strength = strength;
		this.controlStrength = controlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.shape = shape;
		this.shapeZ = shapeZ;
		this.mixProb = mixProb;
		if (mixProb==1)
			unaryEvent = 1;
		this.nearestGene = nearestGene;
		this.distance = distance;
	}
	public GPSPeak(Genome g, String chr, int pos, double strength, 
			double controlStrength, double qvalue, double pvalue, double shape, 
			int unaryEvent, String nearestGene, int distance){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = strength;
		this.controlStrength = controlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.shape = shape;
		this.unaryEvent = unaryEvent;
		if (unaryEvent==1)
			mixProb = 1;
		this.nearestGene = nearestGene;
		this.distance = distance;
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
		return shape;
	}
	public double getShapeZ() {
		return shapeZ;
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
		this.shape = shape;
	}
	public void setShapeZ(double shapeZ) {
		this.shapeZ = shapeZ;
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
	public double getMixProb() {
		return mixProb;
	}
	public void setMixProb(double mixProb) {
		this.mixProb = mixProb;
	}
	public int compareByPValue(GPSPeak p) {
		double diff = getPvalue()- p.getPvalue();
		return	diff==0?0:(diff<0)?1:-1;	//p-value: descending
	}	
	public String toGPS(){
		return toString()+"\t"+"\t"+strength+"\t"+controlStrength+"\t"+qvalue+"\t"+pvalue
		+"\t"+shape+"\t"+unaryEvent+"\t"+nearestGene+"\t"+distance;
	}
	
	public static String toGPS_Header(){
		return "GPS Peak\tIP reads\tControl reads\tQ-value(-log10)\tP-value(-log10)\t"+
		"shape\tunaryEvent\tNearestGene\tDistance";
	}
	
	public String toGPS_short(){
		return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
		+controlStrength+"\t"+qvalue+"\t"+pvalue+"\t"+shape+"\t"+unaryEvent;
	}
	
	public static String toGPS_short_Header(){
		return "Event Location\tNearestGene\tDistance\tIP reads\tControl reads\tQ-value(-log10)\t"+
		"P-value(-log10)\tshape\tunaryEvent";
	}
	
	public String toGPS_motifShort(){
	  return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
	  +controlStrength;//+"\t"+qvalue+"\t"+shape+"\t"+shapeZ;
	}
}
