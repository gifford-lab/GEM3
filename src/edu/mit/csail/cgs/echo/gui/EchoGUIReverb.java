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
import edu.mit.csail.cgs.ewok.types.*;

public class EchoGUIReverb extends EchoComponent {

    private Reverb reverb;
    private Vector<EchoEdge> cachedEdges;
    
    public EchoGUIReverb(EchoGUI g, String n, Rectangle r, Reverb o) {
        super(n, r, g);
        reverb = o;
        reverb.setPeer(this);
        cachedEdges = new Vector<EchoEdge>();
        
        JMenuItem item;
        
        if(reverb.getVerb() instanceof DefaultConstantsParameterized) { 
            popupMenu.add(item = new JMenuItem("Add Default Constants"));
            item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    addDefaultConstants();
                }
            });        
        }
       
        popupMenu.add(item = new JMenuItem("Show Inputs"));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String[] names = reverb.getVerb().getInputNames();
                EchoType[] classes = reverb.getVerb().getInputClasses();
                JDialog dlg = new InformationDialog(new ParameterDisplayPanel(names, classes));
            }
        });        

        popupMenu.add(item = new JMenuItem("Show Output"));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	String[] names;
            	EchoType[] classes;
            	SelfDescribingVerb verb = reverb.getVerb();
        		names = new String[1];
        		classes = new EchoType[1];

            	if(verb.getOutputClass() != null) {
            		names[0] = "Output";
            		classes[0] = reverb.getVerb().getOutputClass();
            	} else { 
            		names[0] = "Output";
            		classes[0] = null;
            	}
            	JDialog dlg = 
            		new InformationDialog(new ParameterDisplayPanel(names, classes));
            }
        });
        
        popupMenu.add(item = new JMenuItem("Show Parameters"));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String[] names = reverb.getVerb().getParameterNames();
                EchoType[] classes = reverb.getVerb().getParameterClasses();
                JDialog dlg = new InformationDialog(new ParameterDisplayPanel(names, classes));
            }
        });
        
    }
    
    public void addDefaultConstants() { 
        SelfDescribingVerb verb = reverb.getVerb();
        if(verb instanceof DefaultConstantsParameterized) { 
            DefaultConstantsParameterized dcp = (DefaultConstantsParameterized)verb;
            String[] knames = dcp.defaultConstantNames();
            SelfDescribingConstant[] kvalues = dcp.defaultConstants();

            EchoGUI gui = getGUI();
            
            for(int i = 0; i < knames.length; i++) { 
                String name = knames[i];
                SelfDescribingConstant k = kvalues[i];
                EchoConstant inp = new EchoConstant(k);
                
                Rectangle rect = gui.getRandomRect(30, 30);
                EchoGUIConstant c = new EchoGUIConstant(gui, k.toString(), rect, inp);
                gui.addEchoComponent(c);
                
                try {
                    reverb.setParam(name, inp);
                } catch (EchoException e) {
                    e.printStackTrace();
                }
            }
            
            rebuildEdges();
        }
    }
    
    public Reverb getEchoObject() { return reverb; }
    
    public void paint(Graphics2D g) { 
        Font oldFont = g.getFont();
        Color oldColor = g.getColor();
        
        int thickness = 3;
        Stroke oldStroke = g.getStroke();
        Stroke thick = new BasicStroke((float)thickness);
        
        Color base = Color.magenta;
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
        g.fillRoundRect(getX(), getY(), getWidth(), getHeight(), curve, curve);
        g.setColor(Color.black);
        g.drawRoundRect(getX(), getY(), getWidth(), getHeight(), curve, curve);
        
        g.setStroke(oldStroke);

        FontRenderContext frc = g.getFontRenderContext();
        TextLayout layout = new TextLayout(getName(), font, frc);
        int textHeight = (int)Math.ceil(layout.getAscent() + layout.getLeading());
        
        g.setColor(Color.black);
        g.drawString(getName(), getX(), getY() + getHeight() + textHeight + 1);
        
        g.setColor(oldColor);
        g.setFont(oldFont);
    }

    public void clearEdges() {
        try {
            SelfDescribingVerb verb = reverb.getVerb();
            
            String[] params = verb.getParameterNames();
            for(int i = 0; params != null && i < params.length; i++) { 
                reverb.clearParam(params[i]);
            }
            
            String[] inputs = verb.getInputNames(); 
            for(int i = 0; inputs != null && i < inputs.length; i++) { 
                reverb.clearInput(inputs[i]);
            }
            
            getGUI().clearOutgoingEdges(this);
            
        } catch (EchoException e) {
            e.printStackTrace();
        } 
        
        rebuildEdges();
    }
    
    public void clearEdge(EchoComponent c) { 
    	if(c instanceof EchoGUIConstant) { 
    		EchoConstant k = ((EchoGUIConstant)c).getEchoObject();
    		reverb.clearParam(k);
    	}
    	
    	if(c instanceof EchoGUIReverb) { 
    		Reverb r = ((EchoGUIReverb)c).getEchoObject();
    		reverb.clearInput(r);
    	}
    	
    	rebuildEdges();
    }

    public void connect(EchoComponent comp) {
        boolean rebuildEdges = false;
        
        if(comp instanceof EchoGUIReverb) { 
            Reverb k = ((EchoGUIReverb)comp).getEchoObject();
            String[] admissibleNames = k.findMatchingInputNames(reverb);
            
            if(admissibleNames != null) { 
                if(admissibleNames.length == 1) { 
                    try {
                        reverb.setInput(admissibleNames[0], k);
                        System.out.println("Matching " + 
                                k.getVerb().getClass().getName() + " --> " + admissibleNames[0]);
                        rebuildEdges = true;
                    } catch (EchoException e) {
                        e.printStackTrace();
                    }
                } else { 
                    AddInputAction action = new AddInputAction(k);
                    SelectionDialog dlg = 
                        new SelectionDialog<String>(new StringArraySelectionPanel(admissibleNames), action);
                    dlg.openDialog();
                }
            }                        
        }
        
        if(comp instanceof EchoGUIConstant) { 
            EchoConstant k = ((EchoGUIConstant)comp).getEchoObject();
            String[] admissibleParams = k.findMatchingParamNames(reverb);
            
            if(admissibleParams != null) { 
                if(admissibleParams.length == 1) { 
                    try {
                        reverb.setParam(admissibleParams[0], k);
                        System.out.println("Matching " + 
                                k.getConstantClass().getName() + " --> " + admissibleParams[0]);
                        rebuildEdges = true;
                    } catch (EchoException e) {
                        e.printStackTrace();
                    }
                } else { 
                    AddParameterAction action = new AddParameterAction(k);
                    SelectionDialog dlg = new SelectionDialog<String>(new StringArraySelectionPanel(admissibleParams), action);
                    dlg.openDialog();
                }
            }
        }
        
        if(rebuildEdges) { rebuildEdges(); } 
    }
    
    public Collection<EchoEdge> getEdges() { return cachedEdges; }
    
    private void rebuildEdges() { 
        cachedEdges = buildEdges();
        repaintGUI();
    }

    private Vector<EchoEdge> buildEdges() {
        Vector<EchoEdge> edges = new Vector<EchoEdge>();
        
        String[] inputNames = reverb.getVerb().getInputNames();
        for(int i = 0; inputNames != null && i < inputNames.length; i++) { 
            Reverb r = reverb.getInput(inputNames[i]);
            if(r != null) { 
                EchoGUIReverb gui = r.getPeer();
                EchoEdge edge = new EchoEdge(gui, this);
                edges.add(edge);
            }
        }

        String[] paramNames = reverb.getVerb().getParameterNames();
        for(int i = 0; paramNames != null && i < paramNames.length; i++) { 
            EchoConstant k = reverb.getParam(paramNames[i]);
            if(k != null) { 
                EchoGUIConstant guik = k.getPeer();
                EchoEdge edge = new EchoEdge(guik, this, paramNames[i]);
                edges.add(edge);
            }
        }
        
        return edges;
    }
    
    private class AddParameterAction implements EchoAction<String> {

        private EchoConstant constant;
        
        public AddParameterAction(EchoConstant k) { 
            constant = k;
        }

        public void doAction(String name) {  
            try {
                reverb.setParam(name, constant);
                rebuildEdges();  
            } catch (EchoException e) {
                e.printStackTrace();
            }
        }
    }

    private class AddInputAction implements EchoAction<String> {

        private Reverb input;
        
        public AddInputAction(Reverb k) { 
            input = k;
        }

        public void doAction(String name) {  
            try {
                reverb.setInput(name, input);
                rebuildEdges();  
            } catch (EchoException e) {
                e.printStackTrace();
            }
        }
    }
}
