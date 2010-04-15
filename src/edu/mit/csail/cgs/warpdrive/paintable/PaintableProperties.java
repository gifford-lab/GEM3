package edu.mit.csail.cgs.warpdrive.paintable;

import edu.mit.csail.cgs.warpdrive.WarpProperties;

public class PaintableProperties extends WarpProperties {

    public Boolean DrawTrackLabel = Boolean.TRUE;
    public String TrackLabel = "";
    
    public PaintableProperties() {
        super();
        DrawTrackLabel = Boolean.TRUE;
        TrackLabel = "";
    }

    public PaintableProperties (boolean drawLabel) {
        super();
        DrawTrackLabel = drawLabel;
        TrackLabel = "";
    }

    public void loadDefaults () {
        // don't load the track label from the defaults, since it varies by experiment.
        String origTrackLabel = TrackLabel;
        super.loadDefaults();
        TrackLabel = origTrackLabel;
    }
    
    public String fileSuffix() {return "wdpp";}
}