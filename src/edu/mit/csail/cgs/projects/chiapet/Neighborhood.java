package edu.mit.csail.cgs.projects.chiapet;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.ProfilePoint;
import edu.mit.csail.cgs.datasets.general.ProfileRegion;
import edu.mit.csail.cgs.datasets.general.Region;

public class Neighborhood {

	private SortedSet<ProfileRegion> profiles;
	
	public float[] getProfile() {
		if (profiles.size()>1) {
			System.err.println("too many profiles");
		}
		ProfileRegion tmp = profiles.first();
		float[] tor = new float[tmp.getEnd()-tmp.getStart()+1];
		for (int i=tmp.getStart(); i<=tmp.getEnd(); i++) {
			tor[i-tmp.getStart()] = (float)tmp.getValue(i);
		}
		return tor;
	}

	public Neighborhood() {
		this.profiles = new TreeSet<ProfileRegion>();
	}

	public void addProfile(ProfileRegion profile) {
		SortedSet<ProfileRegion> headset = profiles.headSet(profile);
		if (headset.size()>0) {
			ProfileRegion last = headset.last();
			if (profile.overlaps(last)) {
				profiles.remove(last);
				profile = profile.combine(last);
			}
		}
		SortedSet<ProfileRegion> tailset = profiles.tailSet(profile);
		if (tailset.size()>0) {
			ProfileRegion first = tailset.first();
			if (profile.overlaps(first)) {
				profiles.remove(first);
				profile = profile.combine(first);
			}
		}
		profiles.add(profile);
	}

	public void collapse() {
		ProfileRegion r = profiles.first();
		for (ProfileRegion pr : profiles) {
			r = r.combine(pr);
		}
		profiles = new TreeSet<ProfileRegion>();
		profiles.add(r);
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		for (ProfileRegion pr : profiles) {
			sb.append(pr.toString()+"\t");
		}
		return sb.toString();
	}

	public String collapsedString() {
		this.collapse();
		ProfileRegion pr = profiles.first();
		StringBuffer sb = new StringBuffer();
		//sb.append(pr.getLocationString());
		for (int i=pr.getStart(); i<pr.getEnd(); i++) {
			sb.append(pr.getValue(i)+"\t");
		}
		return sb.toString();
	}

	public ProfilePoint maxPoint() {
		ProfilePoint max = new ProfilePoint(null, "-1", -1, 0);
		for (ProfileRegion pr : profiles) {
			for (int i=pr.getStart(); i<=pr.getEnd(); i++) {
				if (pr.getValue(i) > max.getProfile()) {
					max = new ProfilePoint(pr.getGenome(), pr.getChrom(), i, pr.getValue(i));
				}
			}
		}
		return max;
	}

	public double getValue(Point p) {
		double tor = 0;
		ProfileRegion r = new ProfileRegion(p.getGenome(), p.getChrom(), p.getLocation(), p.getLocation()+1, new double[2]);
		SortedSet<ProfileRegion> headset = profiles.headSet(r);
		if (headset.size()>0) {
			ProfileRegion last = headset.last();
			if (last.contains(p)) {
				tor = last.getValue(p.getLocation());
			}
		}
		SortedSet<ProfileRegion> tailset = profiles.tailSet(r);
		if (tailset.size()>0) {
			ProfileRegion first = tailset.first();
			if (first.contains(p)) {
				if (tor>0) {
					System.err.println("point contained multiple times in neighborhood");
				}
				tor = first.getValue(p.getLocation());
			}
		}
		return tor;
	}

	public double getTotal() {
		double tor = 0;
		for (ProfileRegion pr : profiles) {
			for (int i=pr.getStart(); i<=pr.getEnd(); i++) {
				tor += pr.getValue(i);
			}
		}
		return tor;
	}
	
	public double getAbsTotal() {
		double tor = 0;
		for (ProfileRegion pr : profiles) {
			tor += pr.getAbsTotal();
		}
		return tor;
	}

	//incomplete
	public double minus(Neighborhood other) {
		double tor = 0;
		boolean bdone = false;
		Iterator<ProfileRegion> itera = this.profiles.iterator();
		Iterator<ProfileRegion> iterb = this.profiles.iterator();
		ProfileRegion a = itera.next();
		ProfileRegion b = iterb.next();
		while (b.getChrom().compareTo(a.getChrom())<0 || (b.getChrom().equals(a.getChrom()) && b.before(a))) {
			for (int i=b.getStart(); i<=b.getEnd(); i++) {
				tor -= b.getValue(i);
			}
			if (!iterb.hasNext()) {
				bdone = true;
				break;
			}
			b = iterb.next();
		}
		if (!iterb.hasNext()) {
			for (int i=a.getStart(); i<=a.getEnd(); i++) {
				tor += a.getValue(i);
			}
		}
		if (b.overlaps(a)) {

		}
		return tor;
	}

	public Neighborhood change(Neighborhood other) {
		Neighborhood tor = new Neighborhood();
		for (ProfileRegion pr : this.profiles) {
			double[] arr = new double[pr.getEnd()-pr.getStart()+1];
			for (int i=pr.getStart(); i<pr.getEnd(); i++) {
				arr[i-pr.getStart()] = pr.getValue(i);
			}
			tor.addProfile(new ProfileRegion(pr.getGenome(), pr.getChrom(), pr.getStart(), pr.getEnd(), arr));
		}
		for (ProfileRegion pr : other.profiles) {
			double[] arr = new double[pr.getEnd()-pr.getStart()+1];
			for (int i=pr.getStart(); i<pr.getEnd(); i++) {
				arr[i-pr.getStart()] = -pr.getValue(i);
			}
			tor.addProfile(new ProfileRegion(pr.getGenome(), pr.getChrom(), pr.getStart(), pr.getEnd(), arr));
		}
		return tor;
	}
	
	public void scale(double f) {
		for (ProfileRegion pr : profiles) {
			pr.scale(f);
		}
	}
	
	public SortedSet<ProfileRegion> topRegions(int binsize, double proportion) {
		SortedSet<ProfileRegion> tor = new TreeSet<ProfileRegion>(new Comparator<ProfileRegion>() {

			public int compare(ProfileRegion arg0, ProfileRegion arg1) {
				double total0 = arg0.getAbsTotal();
				double total1 = arg1.getAbsTotal();
				if (total0>total1) {
					return 1;
				} else if (total0<total1) {
					return -1;
				} else {
					return 0;
				}
			}
			
		});
		double total = this.getAbsTotal();
		double goal = proportion*total;
		double tortal = 0;
		System.err.println(total);
		for (ProfileRegion pr : profiles) {
			for (int i=pr.getStart(); i<pr.getEnd(); i += binsize) {
				ProfileRegion tmp = pr.subRegion(i, i+binsize-1);
				if (tortal<goal) {
					tor.add(tmp);
					tortal += tmp.getAbsTotal();
				} else if (tmp.getAbsTotal()>tor.first().getAbsTotal()) {
					while (tortal>goal) {
						ProfileRegion first = tor.first();
						tor.remove(first);
						tortal -= first.getAbsTotal();
					}
					tor.add(tmp);
					tortal += tmp.getAbsTotal();
				}
			}
		}
		return tor;
	}

}
