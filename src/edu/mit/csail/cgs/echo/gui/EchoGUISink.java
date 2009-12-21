/*
 * Created on Feb 19, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.util.*;

import edu.mit.csail.cgs.echo.*;

public class EchoGUISink extends EchoGUIReverb {
    
    private static final double cos_pi4 = Math.cos(Math.PI / 4.0);
    private static final double pi4 = Math.PI / 4.0;

    public EchoGUISink(EchoGUI g, String n, Rectangle r, Reverb o) {
        super(g, n, r, o);
    }
    
    public void paint(Graphics2D g) { 
        Font oldFont = g.getFont();
        Color oldColor = g.getColor();
        Stroke oldStroke = g.getStroke();
        Stroke thick = new BasicStroke((float)3.0);
        
        Color base = Color.red;
        Color color = params.containsKey("Color") ? (Color)params.get("Color") : base;

        if(isSelected()) { 
            color = params.containsKey("SelectedColor") ?
                    (Color)params.get("SelectedColor") : EchoComponent.lighten(base);
        }
        
        Point cp = getCenterPoint();
        
        Font font = params.containsKey("Font") ? (Font)params.get("Font") : 
            new Font("Times New Roman", Font.PLAIN, 18);
        
        g.setFont(font);

        g.translate(cp.x, cp.y);
        g.rotate(-pi4);

        int nw = (int)Math.round(cos_pi4 * (double)getWidth());
        int nh = (int)Math.round(cos_pi4 * (double)getHeight());
        
        g.setStroke(thick);
        
        g.setColor(color);        
        g.fillRect(-nw/2, -nh/2, nw, nh);
        g.setColor(Color.black);        
        g.drawRect(-nw/2, -nh/2, nw, nh);
        
        g.setStroke(oldStroke);
        
        g.rotate(pi4);
        g.translate(-cp.x, -cp.y);

        FontRenderContext frc = g.getFontRenderContext();
        TextLayout layout = new TextLayout(getName(), font, frc);
        int textHeight = (int)Math.ceil(layout.getAscent() + layout.getLeading());
        
        g.setColor(Color.black);
        g.drawString(getName(), getX(), getY() + getHeight() + textHeight + 1);
        
        g.setColor(oldColor);
        g.setFont(oldFont);
    }
}
