package edu.mit.csail.cgs.viz.paintable;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.awt.Dimension;
import java.awt.Color;
import java.io.File;
import java.io.Writer;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFileChooser;


import edu.mit.csail.cgs.viz.components.ImageConfigurationFrame;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.dom.GenericDOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.DOMImplementation;


public abstract class AbstractPaintable 
implements Paintable, PaintableChangedListener { 
    
	private boolean imageRaster=true, disableEventPassthrough=false;
    public static int sImageWidth, sImageHeight;

    static { 
        sImageWidth = 1000;
        sImageHeight = 750;
    }
    
    protected LinkedList<PaintableChangedListener> fListeners;
    protected int siw, sih;
    
    public AbstractPaintable() {
        fListeners = new LinkedList<PaintableChangedListener>();
        sih = sImageHeight;
        siw = sImageWidth;
    }

    public AbstractPaintable(int h, int w) {
        fListeners = new LinkedList<PaintableChangedListener>();
        sih = h;
        siw = w;
    }
    
    public boolean getImageRaster(){return(imageRaster);}
    public void setImageRaster(boolean ir){imageRaster=ir;}
    public void setImageWidth(int iw) { siw = iw; }
    public void setImageHeight(int ih) { sih = ih; }

    public abstract void paintItem(Graphics g, int x1, int y1, int x2, int y2);
    
    public Collection<Action> getPaintableActions() { 
        LinkedList<Action> lst = new LinkedList<Action>();
        lst.addLast(getSaveImageAction());
        lst.addLast(getImageSettingsAction());
        return lst;
    }
    
    public Action getSaveImageAction() {
    	return(getSaveImageAction("Save Image"));
    }

    public Action getSaveImageAction(String message) {
        return new AbstractAction(message) { 
            /**
             * Comment for <code>serialVersionUID</code>
             */
            private static final long serialVersionUID = 1L;

            public void actionPerformed(ActionEvent e) { 
                String pwdName = System.getProperty("user.dir");
                JFileChooser chooser;
                if(pwdName != null) { 
                    chooser = new JFileChooser(new File(pwdName));
                } else {
                    chooser = new JFileChooser();
                }
                
                int v = 
                    chooser.showSaveDialog(null);
                if(v == JFileChooser.APPROVE_OPTION) { 
                    File f = chooser.getSelectedFile();
                    try {
                        saveImage(f, sImageWidth, sImageHeight, imageRaster);
                        //System.out.println("Saved Image [" + sImageWidth + " by " + sImageHeight +  "]");
                    } catch(IOException ie) {
                        ie.printStackTrace(System.err);
                    }
                }
                
            }
        };
    }
    
    public Action getImageSettingsAction(){
    	final AbstractPaintable p = this;
    	return new AbstractAction("Image Settings..."){
    		public void actionPerformed(ActionEvent e) {
    			new ImageConfigurationFrame(p);
    		}
    	};
    }	
    
    
    public void registerClick(double xf, double yf) { 
    }
    
    protected void setEventPassthrough(boolean v) { 
    	disableEventPassthrough = !v;
    }
    
    public void paintableChanged(PaintableChangedEvent evt) {
    	if(!disableEventPassthrough) { 
    		dispatchChangedEvent();
    	}
    }
    
    public void addPaintableChangedListener(PaintableChangedListener l) { 
        fListeners.addLast(l);
    }
    
    public void removePaintableChangedListener(PaintableChangedListener l) { 
        fListeners.remove(l);
    }
    
    public Image createImage(int w, int h) { 
        BufferedImage im = 
            new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        paintItem(g, 0, 0, w, h);
        return im;
    }
    
    public void saveImage(File f, int w, int h, boolean raster) 
    throws IOException { 

        if (raster) {
            BufferedImage im = 
                new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
            Graphics g = im.getGraphics();
            Graphics2D g2 = (Graphics2D)g;
            g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
            paintItem(g, 0, 0, w, h);
            ImageIO.write(im, "png", f);
        } else {
            DOMImplementation domImpl =
                GenericDOMImplementation.getDOMImplementation();
            // Create an instance of org.w3c.dom.Document
            Document document = domImpl.createDocument(null, "svg", null);
            // Create an instance of the SVG Generator
            SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
            svgGenerator.setSVGCanvasSize(new Dimension(w,h));
            // Ask the test to render into the SVG Graphics2D implementation
            svgGenerator.setColor(Color.white);        
            svgGenerator.fillRect(0,0,w,h);
            paintItem(svgGenerator,25,25,w-50,h-50);

            // Finally, stream out SVG to the standard output using UTF-8
            // character to byte encoding
            boolean useCSS = true; // we want to use CSS style attribute
            Writer out = new OutputStreamWriter(new FileOutputStream(f), "UTF-8");
            svgGenerator.stream(out, useCSS);
        }
    }
    
    protected void dispatchChangedEvent() { 
        PaintableChangedEvent evt = new PaintableChangedEvent(this);
        for(PaintableChangedListener l : fListeners) { 
            l.paintableChanged(evt);
        }
    }    
}
