package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class ProfileRegion extends Region {

	private double[] profile;
	private double total;
	private double abstotal;

	public ProfileRegion(Genome g, String c, int start, int end, double[] profile) {
		super(g, c, start, end);
		this.profile = profile;
		this.total = this.computeTotal();
		this.abstotal = this.computeAbsTotal();
	}

	public ProfileRegion combine(ProfileRegion o) {
		if (!getChrom().equals(o.getChrom())) {
			throw new IllegalArgumentException(getChrom() + " != " + o.getChrom());
		}
		int ns = Math.min(getStart(), o.getStart());
		int ne = Math.max(getEnd(), o.getEnd());
		double[] np = new double[ne-ns+1];
		for (int i=ns; i<=ne; i++) {
			np[i-ns] = this.getValue(i) + o.getValue(i);
		}
		return new ProfileRegion(getGenome(), getChrom(), ns, ne, np);
	}
	
	public double getValue(int position) {
		if (position>=this.getStart() && position<=this.getEnd()) {
			return profile[position-this.getStart()];
		} else {
			return 0;
		}
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(this.getLocationString());
		for (int i=0; i<profile.length; i++) {
			sb.append("|"+profile[i]);
		}
		return sb.toString();
	}
	
	public static ProfileRegion fromString(Genome genome, String s) {
		String[] split = s.split("|");
		Region reg = Region.fromString(genome, split[0]);
		double[] profile = new double[split.length-1];
		for (int i=1; i<profile.length; i++) {
			profile[i-1] = Double.parseDouble(split[i]);
		}
		return new ProfileRegion(genome, reg.getChrom(), reg.getStart(), reg.getEnd(), profile);
	}
	
	public double minus(ProfileRegion o) {
		double tor = 0;
		if (!getChrom().equals(o.getChrom())) {
			throw new IllegalArgumentException(getChrom() + " != " + o.getChrom());
		}
		int ns = Math.min(getStart(), o.getStart());
		int ne = Math.max(getEnd(), o.getEnd());
		for (int i=ns; i<=ne; i++) {
			tor += this.getValue(i) - o.getValue(i);
		}
		return tor;
	}
	
	public ProfileRegion diff(ProfileRegion o) {
		if (!getChrom().equals(o.getChrom())) {
			throw new IllegalArgumentException(getChrom() + " != " + o.getChrom());
		}
		int ns = Math.min(getStart(), o.getStart());
		int ne = Math.max(getEnd(), o.getEnd());
		double[] np = new double[ne-ns+1];
		for (int i=ns; i<=ne; i++) {
			np[i-ns] = this.getValue(i) - o.getValue(i);
		}
		return new ProfileRegion(getGenome(), getChrom(), ns, ne, np);
	}
	
	public void scale(double f) {
		for (int i=0; i<profile.length; i++) {
			profile[i] = f*profile[i];
		}
	}
	
	private double computeTotal() {
		double tor = 0;
		for (int i=0; i<profile.length; i++) {
			tor += profile[i];
		}
		return tor;
	}
	
	public double getTotal() {
		return total;
	}
	
	private double computeAbsTotal() {
		double tor = 0;
		for (int i=0; i<profile.length; i++) {
			tor += Math.abs(profile[i]);
		}
		return tor;
	}
	
	public double getAbsTotal() {
		return abstotal;
	}
	
	public ProfileRegion subRegion(int start, int end) {
		double[] np = new double[end-start+1];
		for (int i=start; i<=end; i++) {
			np[i-start] = this.getValue(i);
		}
		return new ProfileRegion(this.getGenome(), this.getChrom(), start, end, np);
	}

}
