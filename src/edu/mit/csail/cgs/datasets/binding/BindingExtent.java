package edu.mit.csail.cgs.datasets.binding;

import edu.mit.csail.cgs.datasets.species.Genome;

/* BindingExtent extends a BindingEvent by representing a larger
   extent.  For example, the BindingRegion may describe just the
   single point of binding, but the Extent represents the larger region
   around that single point that was incorporated into.  For example, the Extent may be
   a region of elevated probability */
public class BindingExtent extends BindingEvent {

    private int estart, eend;
    
    public BindingExtent(Genome g, String c, int start, int end,
                         double size, double conf, String type,
                         int extentstart, int extentend) {
        super(g,c,start,end,size,conf,type);
        estart = extentstart;
        eend = extentend;
    }
    
    public BindingExtent(BindingEvent evt, int es, int ee) { 
        super(evt);
        estart = es;
        eend = ee;
    }
    
    public BindingExtent(BindingEvent evt) { 
        super(evt);
        estart = evt.getStart();
        eend = evt.getEnd();
    }

    public int getExtentStart() {return estart;}
    public int getExtentEnd() {return eend;}
    
    public boolean equals(Object o) { 
        if(!(o instanceof BindingExtent)) { return false; }
        BindingExtent ext = (BindingExtent)o;
        if(!super.equals(ext)) { return false; }
        if(estart != ext.estart) { return false; }
        if(eend != ext.eend) { return false; }
        return true;
    }

    public int hashCode() { 
        int code = super.hashCode(); 
        code += estart; code *= 37;
        code += eend; code *= 37;
        return code;
    }
}
