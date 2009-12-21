package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.Named;
import edu.mit.csail.cgs.datasets.general.Typed;
import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

public class NamedTypedPainter extends NamedStrandedPainter {
    private RegionExpanderModel model;    
    private Color palegreen = new Color((float)0.0,(float)1.0,(float)0.0,(float).8);
    public NamedTypedPainter(RegionExpanderModel model) {
        super(model);
        this.model = model;
    }
    public Color getColor(Region r) {
        if (r instanceof Typed) {
            String type = ((Typed)r).getType();
            if (type.matches(".*Dubious.*")) {
                return Color.PINK;
            }else if (type.equals("Transposon") || type.equals("repeat_region")) {
                return Color.BLUE;
            } else if (type.matches(".*[rR][nN][aA].*")) {
                return Color.RED;
            } else if (type.equals("Telomeric Region")) {
                return Color.CYAN;
            } else if (type.matches(".*Repeat.*")) {
                return Color.GRAY;
            } else {
                return palegreen;
            }
        } else {
            return super.getColor(r);
        }
    }
    
    public String getLabel(Region r) {
        if (r instanceof Named &&
            r instanceof Typed) {
            return ((Named)r).getName() + "/" + 
                ((Typed)r).getType();
        } else {
            return super.getLabel(r);
        }
    }
}
