package edu.mit.csail.cgs.ewok.verbs.chipseq;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class GPSPeak extends Point{
	double strength;
	double controlStrength;
	double qvalue;
	double pvalue;
	double shape;
	private boolean jointEvent;		// 1 for joint event, 0 for unary, etc
	String nearestGene;
	int distance;
	Point EM_position;

	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double controlStrength, double qvalue, double pvalue, double shape, 
			int joint, String nearestGene, int distance){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = controlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.shape = shape;
		this.jointEvent = joint==1;
		this.nearestGene = nearestGene;
		this.distance = distance;
	}
	
	public GPSPeak(Genome g, String chr, int pos, double ipStrength, 
			double ctrlStrength, double qvalue, double pvalue, double shape){
		super(g, chr.replaceFirst("chr", ""), pos);
		this.strength = ipStrength;
		this.controlStrength = ctrlStrength;
		this.qvalue = qvalue;
		this.pvalue = pvalue;
		this.shape = shape;
	}
	
	 public GPSPeak(Genome g, String chr, int pos, int EM_pos, double strength, 
      double controlStrength, double qvalue, double shape, double shapeZ,
      double mixProb, String nearestGene, int distance){
    super(g, chr.replaceFirst("chr", ""), pos);
    EM_position = new Point(g, chr.replaceFirst("chr", ""), EM_pos);
    this.strength = strength;
    this.controlStrength = controlStrength;
    this.qvalue = qvalue;
    this.shape = shape;
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
		if (mixProb==1)
			jointEvent = false;
		this.nearestGene = nearestGene;
		this.distance = distance;
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
		return shape;
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
		return	diff==0?0:(diff<0)?1:-1;	//p-value: descending
	}
	public int compareByIPStrength(GPSPeak p) {
		double diff = getStrength()- p.getStrength();
		return	diff==0?0:(diff<0)?1:-1;	//ip strength: descending
	}
	public String toGPS(){
		return toString()+"\t"+strength+"\t"+controlStrength+"\t"+qvalue+"\t"+pvalue
		+"\t"+shape+"\t"+jointEvent+"\t"+nearestGene+"\t"+distance;
	}
	
	public static String toGPS_Header(){
		return "GPS Event\tIP\tControl\tQ_-lg10\tP_-lg10\t"+
		"Shape\tJoint\tNearestGene\tDistance";
	}
	
	public String toGPS_short(){
		return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
		+controlStrength+"\t"+qvalue+"\t"+pvalue+"\t"+shape+"\t"+jointEvent;
	}
	
	public static String toGPS_short_Header(){
		return "Event Location\tNearestGene\tDist\tIP\tControl\tQ_-lg10\t"+
		"P_-lg10\tShape\tJoint";
	}
	
	public String toGPS_motifShort(){
	  return toString()+"\t"+nearestGene+"\t"+distance+"\t"+strength+"\t"
	  +controlStrength;//+"\t"+qvalue+"\t"+shape+"\t"+shapeZ;
	}
}
