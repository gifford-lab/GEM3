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
import edu.mit.csail.cgs.ewok.types.EchoType;
import edu.mit.csail.cgs.ewok.types.SelfDescribingConstant;
import edu.mit.csail.cgs.ewok.types.ValueWrapper;

public class EchoGUIConstant extends EchoComponent implements SelectionDialogListener {
    
    private static Vector<EchoEdge> emptyEdges = new Vector<EchoEdge>();

    private EchoConstant echoConst;
    private SelectionDialog changeDlg;
    
    public EchoGUIConstant(EchoGUI g, String n, Rectangle r, EchoConstant o) {
        super(n, r, g);
    
        echoConst = o;
        echoConst.setPeer(this);
        
        SelfDescribingConstant k = echoConst.getConstant();
        if(!(k instanceof ValueWrapper)) { throw new IllegalStateException(); }
        ValueWrapper vw = (ValueWrapper)k;
        changeDlg = SelectionDialog.fromValueWrapper(vw);

        if(changeDlg != null) {
            changeDlg.addSelectionDialogListener(this);
            
            JMenuItem item;

            popupMenu.add(item = new JMenuItem("Edit Value"));
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) { 
                    changeDlg.openDialog();
                }
            });
        }
    }
    
    public void selectionDialogClosed(SelectionDialogEvent evt) { 
        if(changeDlg.isAcceptedChoice()) { 
            Object val = changeDlg.getSelectedValue();
            setName(val.toString());
            repaintGUI();
        }
    }
    
    public EchoConstant getEchoObject() { return echoConst; }
    
    public void paint(Graphics2D g) { 
        Font oldFont = g.getFont();
        Color oldColor = g.getColor();
        Stroke oldStroke = g.getStroke();
        Stroke thick = new BasicStroke((float)3.0);
        
        Color base = Color.blue;
        Color color = params.containsKey("Color") ? (Color)params.get("Color") : base;
        if(isSelected()) { 
            color = params.containsKey("SelectedColor") ?
                    (Color)params.get("SelectedColor") : EchoComponent.lighten(base);
        }
        
        Font font = params.containsKey("Font") ? (Font)params.get("Font") : 
            new Font("Times New Roman", Font.PLAIN, 18);
        
        g.setFont(font);

        g.setStroke(thick);
        
        g.setColor(color);
        g.fillOval(getX(), getY(), getWidth(), getHeight());
        g.setColor(Color.black);
        g.drawOval(getX(), getY(), getWidth(), getHeight());
        
        g.setStroke(oldStroke);

        FontRenderContext frc = g.getFontRenderContext();
        TextLayout layout = new TextLayout(getName(), font, frc);
        int textHeight = (int)Math.ceil(layout.getAscent() + layout.getLeading());
        
        g.setColor(Color.black);
        g.drawString(getName(), getX(), getY() + getHeight() + textHeight + 1);
        
        g.setColor(oldColor);
        g.setFont(oldFont);
    }
    
    public boolean couldParameterize(EchoComponent comp) {
    	if(comp instanceof EchoGUIReverb) { 
    		Reverb rev = ((EchoGUIReverb)comp).getEchoObject();
    		EchoType[] pc = rev.getVerb().getParameterClasses();
    		for(int i = 0; pc != null && i < pc.length; i++) {
    			if(echoConst.getConstantClass().isSubType(pc[i])) { 
    				return true;
    			}
    		}    		
    	}
        
    	return false;
    }

    public void clearEdges() {
    }
    
    public void clearEdge(EchoComponent c) { 
    }

    public void connect(EchoComponent comp) {
    	if(!(comp instanceof EchoGUIConstant)) { 
    		comp.connect(this);
    	}
    }

    public Collection<EchoEdge> getEdges() {
        return emptyEdges;
    }
}
