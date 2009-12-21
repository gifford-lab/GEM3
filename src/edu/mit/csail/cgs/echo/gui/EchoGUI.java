/*
 * Created on Feb 19, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;

import javax.swing.JPanel;
import javax.swing.JPopupMenu;

import edu.mit.csail.cgs.echo.EchoBase;
import edu.mit.csail.cgs.echo.EchoConstant;
import edu.mit.csail.cgs.echo.EchoException;
import edu.mit.csail.cgs.echo.Reverb;

public class EchoGUI extends JPanel {
    
    private EchoBase base;
    private Vector<EchoComponent> comps;
    
    private int movingIndex;
    private Point lastMoveLocation;
    
    public EchoGUI() { 
        super();
        comps = new Vector<EchoComponent>();
        base = new EchoBase();

        movingIndex = -1;
        lastMoveLocation = null;
        
        addMouseListener(new MouseAdapter() {

            public void mouseClicked(MouseEvent me) {
                Vector<Integer> cp = getComponentsAtPoint(me.getPoint());

                if(me.getButton() == MouseEvent.BUTTON3 || me.isControlDown()) { 
                	if(cp.size() > 0) { 
                		int idx = cp.get(0);
                		JPopupMenu menu = comps.get(idx).getPopupMenu();
                		if(menu != null) {
                			System.out.println("Showing menu...");
                			menu.show(EchoGUI.this, me.getX(), me.getY()); 
                		}                   
                	}
                }
            }

            public void mousePressed(MouseEvent me) {
                Vector<Integer> cp = getComponentsAtPoint(me.getPoint());
                
                if(me.getButton() == MouseEvent.BUTTON1) { 
                    if(cp.size() > 0) { 
                        movingIndex = cp.get(0);
                        lastMoveLocation = me.getPoint();
                        
                        EchoComponent comp = comps.get(movingIndex);
                        comp.setSelected(true);
                        setParameterizedSelected(comp);
                    }
                }
            }

            public void mouseReleased(MouseEvent me) {
                if(me.getButton() == MouseEvent.BUTTON1) { 
                    if(movingIndex != -1) { 
                        searchAndSetInput(movingIndex);

                        EchoComponent rem = comps.remove(movingIndex);
                        comps.insertElementAt(rem, 0);
                        clearSelected();
                    }
                    movingIndex = -1;
                    lastMoveLocation = null;
                }
            } 
            
        });
        
        addMouseMotionListener(new MouseMotionAdapter() {
            public void mouseDragged(MouseEvent me) {
                if(lastMoveLocation != null) { 
                    int dx = me.getPoint().x - lastMoveLocation.x;
                    int dy = me.getPoint().y - lastMoveLocation.y;
                    comps.get(movingIndex).move(dx, dy);
                    lastMoveLocation = me.getPoint(); 
                    repaint();
                }
            }
        });
    }
    
    public Dimension getPreferredSize() { return new Dimension(400, 400); }
    
    public void clearOutgoingEdges(EchoComponent c) { 
    	for(EchoComponent comp : comps) { 
    		comp.clearEdge(c);
    	}
    }
    
    public void clear() {
        base.clear();
        comps.clear();
        repaint();
    }
    
    public void remove(EchoComponent c) {
    	if(c instanceof EchoGUIReverb) {
    		EchoGUIReverb r = (EchoGUIReverb)c;
    		base.remove(r.getEchoObject());
    	}
    	
    	if(c instanceof EchoGUIConstant) { 
    		EchoGUIConstant k = (EchoGUIConstant)c;
    		base.remove(k.getEchoObject());
    	}
    	
    	c.clearEdges();
    	for(EchoComponent comp : comps) { 
    		comp.clearEdge(c);
    	}
        	
    	comps.remove(c);
    	repaint();
    }

    private synchronized void searchAndSetInput(int cindex) {
        EchoComponent comp1 = comps.get(cindex);
        
        for(int i = 0; i < comps.size(); i++) { 
            if(i != cindex) { 
                EchoComponent comp2 = comps.get(i);
                if(comp1.containsPoint(comp2.getCenterPoint())) { 
                    comp1.connect(comp2);
                    return;
                }
            }
        }
    }
    
    public Rectangle getRandomRect(int w, int h) {
        Random rand = new Random();
        int x = rand.nextInt(getWidth()-w);
        int y = rand.nextInt(getHeight()-h);
        return new Rectangle(x, y, w, h);
    }
    
    private Vector<Integer> getComponentsAtPoint(Point p) { 
        Vector<Integer> v = new Vector<Integer>();
        synchronized(this) { 
            for(int i = 0; i < comps.size(); i++) { 
                if(comps.get(i).containsPoint(p)) { 
                    v.add(i);
                }
            }
        }
        return v;
    }
    
    public synchronized void addEchoComponent(EchoComponent ec) {
        comps.add(ec);
        if(ec instanceof EchoGUIConstant) {
            EchoGUIConstant v = (EchoGUIConstant)ec;
            base.addConstant(v.getEchoObject());
        }
        if(ec instanceof EchoGUIReverb) {
            EchoGUIReverb v = (EchoGUIReverb)ec;
            base.addReverb(v.getEchoObject());
        }
        repaint();
    }
    
    public void runBase() { 
        try {
            base.runEchoProgram();
        } catch (EchoException e) {
            e.printStackTrace(System.err);
        }
    }
    
    protected void paintComponent(Graphics g) { 
        super.paintComponent(g);
        
        Graphics2D g2 = (Graphics2D)g;
        Map oldHints = g2.getRenderingHints();
        
        Map newHints = new HashMap(oldHints);
        newHints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHints(newHints);
        
        int w = getWidth(), h = getHeight();
        g.setColor(Color.white);
        g.fillRect(0, 0, w, h);
        
        g.setColor(Color.black);
        g.drawRect(0, 0, w-1, h-1);
        
        synchronized(this) { 
            for(int i = 0; i < comps.size(); i++) {
                for(EchoEdge ee : comps.get(i).getEdges()) { 
                    ee.paint(g2);
                }
            }
            
            for(int i = comps.size()-1 ; i >= 0; i--) { 
                EchoComponent ec = comps.get(i);
                ec.paint(g2);
            }
        }
        
        g2.setRenderingHints(oldHints);
    }
    
    public void setParameterizedSelected(EchoComponent comp) {
    	if(comp instanceof EchoGUIConstant) {
    		EchoGUIConstant k = (EchoGUIConstant)comp;
    		for(int i = 0; i < comps.size(); i++) { 
    			if(k.couldParameterize(comps.get(i))) { 
    				comps.get(i).setSelected(true);
    			}
    		}
    	} else { 
    		for(int i = 0; i < comps.size(); i++) {
    			if(comps.get(i) instanceof EchoGUIConstant) { 
    				EchoGUIConstant k = (EchoGUIConstant)comps.get(i);
    				if(k.couldParameterize(comp)) { 
    					k.setSelected(true);
    				}
    			}
    		}
    	}
    	repaint();
    }
    
    public void clearSelected() { 
    	for(int i = 0; i < comps.size(); i++) { comps.get(i).setSelected(false); }
    	repaint();
    }

    private class AddEdgeAction implements EchoAction<String> {
    	
    	private EchoComponent source, target;
    	private EchoConstant param;
    	private Reverb inpInter;
    	private EchoEdgeConnector edgeTracker;
    	
    	public AddEdgeAction(EchoComponent s, EchoComponent targ, 
    			EchoEdgeConnector t, EchoConstant k, Reverb eii) { 
    		source = s;
    		target = targ;
    		edgeTracker = t;
    		inpInter = eii;
    		param = k;
    	}
    	
    	public void doAction(String name) { 
    		try {
				inpInter.setParam(name, param);
	    		EchoEdge edge = new EchoEdge(source, target, name);
	    		edgeTracker.connect(source);
	    		EchoGUI.this.repaint();
			} catch (EchoException e) {
				e.printStackTrace();
			}
    	}
    }
}

