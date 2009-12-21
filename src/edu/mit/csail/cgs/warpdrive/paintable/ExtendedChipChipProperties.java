package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.Color;

public class ExtendedChipChipProperties extends ChipChipProperties {

    public Boolean IntensitiesOnLogScale = Boolean.TRUE;
    public Integer AverageIntensitiesAcrossNPoints = 1;
    public Integer AverageRatiosAcrossNPoints = 1;
    public Boolean AverageAcrossReplicates = Boolean.FALSE;
    public Integer MinLineWidth = 2;
    public Integer MinCircleRadius = 2;
    public Integer FontSize = 8;
    public Boolean DrawLabelOnRight = Boolean.TRUE;
    public Boolean DrawCy5 = Boolean.FALSE;
    public Boolean DrawCy3 = Boolean.FALSE;
    public Boolean DrawStrandedRatio = Boolean.FALSE;
    public Double MinIntensity = 10.0;
    public Double MaxIntensity = 70000.0;
    public Integer MaximumLineDistance = 500;
    public Color Color = java.awt.Color.BLUE;
}