package edu.mit.csail.cgs.warpdrive.paintable;

public class ChipSeqAnalysisProperties extends PaintableProperties {

    public void loadDefaults () {
        // don't load the track label from the defaults, since it varies by experiment.
        String origTrackLabel = TrackLabel;
        super.loadDefaults();
        TrackLabel = origTrackLabel;
    }

}