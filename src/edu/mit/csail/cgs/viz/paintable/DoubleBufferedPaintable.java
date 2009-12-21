/*
 * Author: tdanford
 * Date: May 18, 2008
 */
/**
 * 
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFileChooser;

/**
 * @author tdanford
 *
 */
public class DoubleBufferedPaintable implements Paintable, PaintableChangedListener {
	
	private Paintable inner;
	private BufferedImage buffer;
	private int width, height;
	private LinkedList<PaintableChangedListener> listeners;
	private boolean paintBackground;
	
	public DoubleBufferedPaintable(Paintable p) { 
		inner = p;
		width = height = -1;
		buffer = null;
		paintBackground = false;
		listeners = new LinkedList<PaintableChangedListener>();
		inner.addPaintableChangedListener(this);
	}
	
	public void paintBackground(boolean bg) { 
		paintBackground = bg;
		dispatchChangedEvent();
	}
	
	public void invalidateBuffer() { 
		buffer = null;
		width = height = -1;
	}
	
	private void createNewBuffer(int w, int h) {
        BufferedImage im = 
            new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, 
        		RenderingHints.VALUE_ANTIALIAS_ON));
    
        if(paintBackground) { 
        	g.setColor(Color.white);
        	g.fillRect(0, 0, w, h);
        }
        
        inner.paintItem(g, 0, 0, w, h);
        
		width = w; height = h;
        buffer = im;
	}

	public void paintableChanged(PaintableChangedEvent e) {
		invalidateBuffer();
		dispatchChangedEvent();
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#addPaintableChangedListener(edu.mit.csail.cgs.viz.paintable.PaintableChangedListener)
	 */
	public void addPaintableChangedListener(PaintableChangedListener l) {
		listeners.add(l);
	}
	
	protected void dispatchChangedEvent() { 
		PaintableChangedEvent evt = new PaintableChangedEvent(this);
		for(PaintableChangedListener listener : listeners) { 
			listener.paintableChanged(evt);
		}
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#getPaintableActions()
	 */
	public Collection<Action> getPaintableActions() {
		return inner.getPaintableActions();
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#getSaveImageAction()
	 */
    public Action getSaveImageAction() { 
        return new AbstractAction("Save As Image...") { 
            public void actionPerformed(ActionEvent e) { 
                String pwdName = System.getProperty("user.dir");
                JFileChooser chooser;
                if(pwdName != null) { 
                    chooser = new JFileChooser(new File(pwdName));
                } else {
                    chooser = new JFileChooser();
                }
                
                int v = 
                    chooser.showOpenDialog(null);
                if(v == JFileChooser.APPROVE_OPTION) { 
                    File f = chooser.getSelectedFile();
                    int w = width > 0 ? width : 800;
                    int h = height > 0 ? height : 600;
                    try {
                        saveImage(f, w, h, true);
                        System.out.println("Saved Image [" + width + " by " + height +  "]");
                    } catch(IOException ie) {
                        ie.printStackTrace(System.err);
                    }
                }
                
            }
        };
    }
    
	protected void saveImage(File f, int w, int h, boolean b) throws IOException {
		if(buffer == null || w != width || h != height) { 
			createNewBuffer(w, h);
		}
        ImageIO.write(buffer, "png", f);
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#paintItem(java.awt.Graphics, int, int, int, int)
	 */
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2 - x1, h = y2 - y1;
		if(buffer == null || w != width || h != height) { 
			createNewBuffer(w, h);
		}
		g.setColor(Color.white);
		g.fillRect(0, 0, w, h);
		g.drawImage(buffer, 0, 0, w, h, null);
	}
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#registerClick(double, double)
	 */
	public void registerClick(double xf, double yf) {
		inner.registerClick(xf, yf);
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.paintable.Paintable#removePaintableChangedListener(edu.mit.csail.cgs.viz.paintable.PaintableChangedListener)
	 */
	public void removePaintableChangedListener(PaintableChangedListener l) {
		listeners.remove(l);
	}
}
