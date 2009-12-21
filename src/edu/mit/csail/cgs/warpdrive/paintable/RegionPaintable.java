package edu.mit.csail.cgs.warpdrive.paintable;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.*;

/* A RegionPaintable displays some dataset across a genomic region */

public abstract class RegionPaintable extends WarpPaintable {

    private Region region;

    public RegionPaintable() {
    }

    /* How much vertical space can this Paintable use? 
       This is returned in pixels and basically assumes screen
       resolution.  -1 means "give me as much space as possible".
    */
    public int getMaxVertSpace () {return -1;}
    public Region getRegion() {return region;}
    /* Sets the genomic region that this Paintable will display.
       The Paintable should  set canPaint() to false
       until it is ready to display.  */
    public void setRegion(Region r) {
        if (r == null) {
            throw new NullPointerException("RegionPaintable won't take a null region");
        }
        region = r;
        setCanPaint(false);
        setWantsPaint(false);
    }
    /* utility functions to help RegionPainters calculate display coordinate.  These take
       region or value boundaries (eg, chromosome starting and ending coordinate or
       maximum and minimum value) and the area into which the paintable is drawing
       and maps a specific coordinate or value to an X or Y value.*/
    public int getXPos(int pos, int start, int end, int leftx, int rightx) {
        if (pos < start) {return leftx;}
        if (pos > end) {return rightx;}
        return (int)((((float)(pos - start))/((float)(end - start))) * (rightx - leftx) + leftx);
    }
    public int getYPos(double val, double minval, double maxval, int miny, int maxy, boolean logScale) {
        if (logScale) {
            return getYPosLog(val,minval,maxval,miny,maxy);
        }
        if (val < minval) {return maxy;}
        if (val > maxval) {return miny;}
        return maxy - (int)(((val - minval)/(maxval - minval)) * ((double)(maxy - miny)));
    }
    public int getYPosLog(double val, double minval, double maxval, int miny, int maxy) {
        if (minval <= .001) {minval = .001;}
        if (val < minval) {return maxy;}
        if (val > maxval) {return miny;}
        val = Math.log(val);
        minval = Math.log(minval);
        maxval = Math.log(maxval);
        return maxy - (int)(((val - minval)/(maxval - minval)) * ((double)(maxy - miny)));
    }
}
