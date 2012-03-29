package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.projects.chiapet.PairedStorage;

public class InteractionLikelihoodModel extends WarpModel implements
RegionModel, Runnable {

	private TreeMap<Integer,Double> results; //map from left edge of bin to height
	private Region region;
	private boolean newinput;
	private InteractionLikelihoodProperties props;
	private PairedStorage storage;
	private String currentChrom;

	public InteractionLikelihoodModel(Genome genome, PairedStorage storage) throws NumberFormatException, IOException {
		props = new InteractionLikelihoodProperties();
		region = null;
		newinput = false;
		this.storage = storage;
		currentChrom = Region.fromString(genome, props.Anchor).getChrom();
		storage.loadChrom(currentChrom);
	}

	public boolean isReady() {return !newinput;}

	public synchronized void run() {
		while(keepRunning()) {
			try {
				if (!newinput) {
					wait();
				}
			} catch (InterruptedException ex) {

			}
			if (newinput) {
				System.err.println(props);
				Region anchor = Region.fromString(region.getGenome(), props.Anchor);
				if (!anchor.getChrom().equals(currentChrom)) {
					currentChrom = anchor.getChrom();
					try {
						storage.loadChrom(currentChrom);
					} catch (NumberFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				int width = props.BinWidth;
				results = new TreeMap<Integer,Double>();
				if (!anchor.getChrom().equals(region.getChrom())) {
					int tmpx = region.getStart();
					System.err.println(region.getWidth()/width);
					while (tmpx<region.getEnd()) {
						Region tmpreg = new Region(region.getGenome(),region.getChrom(),tmpx,tmpx+width);
						double ratio = Math.max(0, storage.computeRatio(anchor, tmpreg));
						results.put(tmpx, ratio);
						tmpx += width;
					}
				} else if (region.getStart()<anchor.getStart()) {
					int tmpx = anchor.getStart();
					while (tmpx>region.getStart()) {
						tmpx -= width;
					}
					while (tmpx<region.getEnd()) {
						Region tmpreg = new Region(region.getGenome(),region.getChrom(),tmpx,tmpx+width);
						double ratio = Math.max(0, storage.computeRatio(anchor, tmpreg));
						results.put(tmpx, ratio);
						tmpx += width;
					}
				} else {
					int tmpx = anchor.getStart();
					while (tmpx<region.getStart()) {
						tmpx += width;
					}
					tmpx -= width;
					while (tmpx<region.getEnd()) {
						Region tmpreg = new Region(region.getGenome(),region.getChrom(),tmpx,tmpx+width);
						double ratio = Math.max(0, storage.computeRatio(anchor, tmpreg));
						results.put(tmpx, ratio);
						tmpx += width;
					}
				}
			}
			System.err.println("computed ratios");
			newinput = false;
			notifyListeners();
		}
	}

	public void setRegion(Region r) {
		if (newinput == false) {
			if (!r.equals(region)) {
				region = r;
				newinput = true;
			} else {
				notifyListeners();
			}
		}
	}

	public Region getRegion() {return region;}

	public Map<Integer,Double> getResults() {
		return results;
	}

	public InteractionLikelihoodProperties getProperties() {
		return props;
	}

}
