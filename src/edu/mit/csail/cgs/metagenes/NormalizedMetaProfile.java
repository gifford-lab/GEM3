package edu.mit.csail.cgs.metagenes;

public class NormalizedMetaProfile extends MetaProfile{

	public NormalizedMetaProfile(String n, BinningParameters bps) {
		super(n, bps);
	}

	public synchronized void addProfile(Profile p) {
		if(p.isStranded()){
			stranded=true;
		}
		if(p.length() != params.getNumBins()) { 
			throw new IllegalArgumentException(String.format("Profile length %d doesn't" +
					" match bin-length %d", p.length(), params.getNumBins()));
		}
		
		if(isNormalized()) { 
			throw new IllegalArgumentException("Can't add profile to a normalized MetaProfile");
		}
		
		if(profiles.contains(p)) { 
			/*throw new IllegalArgumentException(String.format(
					"Can't add same profile %s to MetaProfile", 
					p.getName()));*/
		}else{
		
			profiles.add(p);
			double count = profiles.size();
			p.addProfileListener(this);
			min=max=0.0;
			for(int i = 0; i< values.length ;i++) { 
				values[i] = (values[i]*(count-1.0)/(count))+(p.value(i)/count);
				max = Math.max(max, values[i]);
				min = Math.min(min, values[i]);
			}
			
			dispatchChange(new ProfileEvent(this, p));
		}
	}

}
