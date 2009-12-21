package edu.mit.csail.cgs.warpdrive.model;

import edu.mit.csail.cgs.datasets.general.Region;

/* A RegionModel is a model that also carries an ewok Region with it.  */
   
public interface RegionModel extends Model {

    public void setRegion(Region r) throws NullPointerException;
    public Region getRegion();
}
