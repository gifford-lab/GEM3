package edu.mit.csail.cgs.warpdrive.model;

import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.projects.readdb.PairedHit;

public class InteractionArcModel extends WarpModel implements RegionModel,
		Runnable {
	
	private List<PairedHit> results, otherchrom;

	public boolean isReady() {
		// TODO Auto-generated method stub
		return false;
	}

	public void run() {
		// TODO Auto-generated method stub

	}

	public void setRegion(Region r) throws NullPointerException {
		// TODO Auto-generated method stub

	}

	public Region getRegion() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public List<PairedHit> getResults () {return results;}

}
