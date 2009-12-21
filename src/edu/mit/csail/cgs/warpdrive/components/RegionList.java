package edu.mit.csail.cgs.warpdrive.components;

import edu.mit.csail.cgs.datasets.general.Region;

public interface RegionList {
    public void addRegion(Region r);
    public int regionListSize();
    public Region regionAt(int i);
}
