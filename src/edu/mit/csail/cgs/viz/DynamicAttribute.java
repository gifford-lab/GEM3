/*
 * Created on Apr 6, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package edu.mit.csail.cgs.viz;

import java.awt.*;

/**
 * @author Lindy Briggs
 * 
 * (Again, another class by Lindy: this was, as with AttributeCalculator, 
 * originally copied over from psrg.lindy package.  This is a default implementation
 * of the AttributeCalculator class).
 * 
 * The attributes in a DynamicAttribute object change depending on the size of the 
 * viewing window, unless the user specifically sets attribute values. If that
 * happens, the specified attributes will remain fixed. They can be reset to their
 * original dynamic version as well. Right now, the dynamic attribute values are 
 * kind of random on my part, for experiment's sake.
 * 
 * I may consider adding a way for the user to set values as percentages later on.
 */
public class DynamicAttribute implements AttributeCalculator {
	
    private static DynamicAttribute global = null;
    public static final int SCREEN = 1,
        PRINT = 2,
        DISPLAY = 3;
    private int ptWidthDiv, lineWidthDiv, fontSizeDiv; 
    private int type = SCREEN;

    public static DynamicAttribute getGlobalAttributes() {
        if (global == null) {
            global = new DynamicAttribute();
        }
        return global;
    }

    public DynamicAttribute() {
        ptWidthDiv = 125;
        lineWidthDiv = 150;
        fontSizeDiv = 75;
    }
    
    public int getType() {return type;}
    public int getPointWidth(int w, int h) { return getPointWidth(w, h, 1.0); }
    public int getLineWidth(int w, int h) { return getLineWidth(w, h, 1.0); }
    public int getFontSize(int w, int h) { return getFontSize(w, h, 1.0); }
    
    public double getTypeScale() {
        if (type == SCREEN) {
            return .8;
        }
        if (type == PRINT) {
            return 1.4;
        }
        if (type == DISPLAY) {
            return 1.2;
        }
        return 1.0;
    }

    // get diameter of a circle (although the agilent painter is treating this as a radius...)
    public int getPointWidth(int w, int h, double scale) {
        int d = Math.min(w, h);        
        int ptWidth = (int)Math.round((d / (double)ptWidthDiv) * scale * getTypeScale());
        if (ptWidth <= 0) { 
            ptWidth = 1;
        }
        return ptWidth;
    }

    // get width of a line
    public int getLineWidth(int w, int h, double scale) {
        int d = Math.min(w, h);
        int lnWidth = (int)Math.round(scale * (d / (double)lineWidthDiv) * getTypeScale());
        if (lnWidth <= 0) { 
            lnWidth = 1;
        }
        return lnWidth;
    }

    // get size of a font
    public int getFontSize(int w, int h, double scale) {
        int d = Math.min(w, h * 15);
        int ftSize = (int)Math.round(scale * ((d / (double)fontSizeDiv) + 10.0));
        if (type == SCREEN) {
            ftSize = Math.max(Math.min(ftSize,18),8);
        }
        if (type == PRINT) {
            ftSize = Math.max(Math.min(ftSize,50),14);
        }
        if (type == DISPLAY) {
            ftSize = Math.max(Math.min(ftSize,36),14);
        }
        return ftSize;
    }
    public Font getLargeLabelFont(int w, int h) {
        return getLargeLabelFont(w,h,1.0);
    }
    public Font getRegionLabelFont(int w, int h) {
        return getRegionLabelFont(w,h,1.0);
    }
    public Font getPointLabelFont(int w, int h) {
        return getPointLabelFont(w,h,1.0);
    }    
    public Font getLargeLabelFont(int w, int h, double scale) {
        int size = getFontSize(Math.max(w,h), Math.max(w,h),1.4 * scale * getTypeScale());
        if (size > Math.min(w,h)) {
            size = Math.min(w,h);
        }
        return new Font("Arial",Font.PLAIN,size);
    }
    public Font getRegionLabelFont(int w, int h, double scale) {
        return new Font("Arial",Font.PLAIN,getFontSize(w,h,1.1 * scale * getTypeScale()));
    }
    public Font getPointLabelFont(int w, int h, double scale) {
        return new Font("Arial",Font.PLAIN,getFontSize(w,h,.8 * scale * getTypeScale()));
    }
    

    public void setFontSize(int size) {
        fontSizeDiv = size;
    }
    public void setLineWidth(int w) {
        lineWidthDiv = w;
    }
    public void setPointWidth(int w) {
        ptWidthDiv = w;
    }
    public void setType(int t) {
        type = t;
    }
	
}
