/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.mit.csail.cgs.metagenes;

public interface Profile {
	public String getName();
	public double value(int i);
	public int length();
	public double max();
	public double min();
	public int getNumProfiles();
	public void setStranded(boolean s);
	public boolean isStranded();
	
	public BinningParameters getBinningParameters();
	
	public void addProfileListener(ProfileListener pl);
	public void removeProfileListener(ProfileListener pl);
}