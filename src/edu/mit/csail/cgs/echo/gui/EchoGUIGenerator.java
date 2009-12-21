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

public class EchoGUIGenerator extends EchoGUIReverb {    
    public EchoGUIGenerator(EchoGUI g, String n, Rectangle r, Reverb o) {
        super(g, n, r, o);
    }
    
    public void paint(Graphics2D g) { 
        Font oldFont = g.getFont();
        Color oldColor = g.getColor();
        
        int thickness = 3;
        Stroke oldStroke = g.getStroke();
        Stroke thick = new BasicStroke((float)thickness);
        
        Color base = Color.green;
        Color color = params.containsKey("Color") ? (Color)params.get("Color") : base; 
        if(isSelected()) { 
            color = params.containsKey("SelectedColor") ?
                    (Color)params.get("SelectedColor") : EchoComponent.lighten(base);
        }
        
        Font font = params.containsKey("Font") ? (Font)params.get("Font") : 
            new Font("Times New Roman", Font.PLAIN, 18);
        
        g.setFont(font);

        g.setStroke(thick);
        
        int curve = thickness*3+2;
        g.setColor(color);
        g.fillRect(getX(), getY(), getWidth(), getHeight());
        g.setColor(Color.black);
        g.drawRect(getX(), getY(), getWidth(), getHeight());
        
        g.setStroke(oldStroke);

        FontRenderContext frc = g.getFontRenderContext();
        TextLayout layout = new TextLayout(getName(), font, frc);
        int textHeight = (int)Math.ceil(layout.getAscent() + layout.getLeading());
        
        g.setColor(Color.black);
        g.drawString(getName(), getX(), getY() + getHeight() + textHeight + 1);
        
        g.setColor(oldColor);
        g.setFont(oldFont);
    }


}
