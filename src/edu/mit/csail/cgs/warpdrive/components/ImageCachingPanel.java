/*
 * Created on Jan 24, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;
import javax.swing.*;

public class ImageCachingPanel extends JPanel {

    private Image cachedImage;
    private JPanel internal;
    
    public ImageCachingPanel(JPanel p) { 
        internal = p;
        setLayout(new BorderLayout());
        add(internal, BorderLayout.CENTER);
        cachedImage = null;
        
        this.addComponentListener(new ComponentAdapter() {
            public void componentResized(ComponentEvent e) {
                invalidateImage();
            }
        });
    }
    
    public void invalidateImage() { 
        cachedImage = null;
        repaint();
    }

    public void paint(Graphics g) { 
        try { 
            int w = getWidth(), h = getHeight();
            if(cachedImage == null) {
                cachedImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
                Graphics ig = cachedImage.getGraphics();
                Graphics2D g2 = (Graphics2D)ig;
                internal.paint(g2);
            }
        
            g.drawImage(cachedImage, 0, 0, w, h, null);
            
        } catch(Exception e) { 
            e.printStackTrace(System.err);
        }
        
        super.paintBorder(g);
    }
}
