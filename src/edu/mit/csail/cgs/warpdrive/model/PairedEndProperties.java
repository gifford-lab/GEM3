package edu.mit.csail.cgs.warpdrive.model;

public class PairedEndProperties extends ModelProperties {

    public Double MinimumDistance=1.0;
    public Boolean DeDuplicateByPosition = false;
    public Boolean LeftAlwaysLesser = true;
    public Boolean PrintData = false;
    public Boolean ShowSelfLigation = true;
    public Boolean RightFlipped = true;
    public Integer SelfLigationCutoff = 10000;

}