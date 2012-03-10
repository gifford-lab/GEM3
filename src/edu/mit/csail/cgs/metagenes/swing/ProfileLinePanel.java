package edu.mit.csail.cgs.metagenes.swing;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;

import javax.imageio.ImageIO;
import javax.swing.*;

import edu.mit.csail.cgs.metagenes.*;
import edu.mit.csail.cgs.viz.paintable.PaintableScale;

public class ProfileLinePanel extends JPanel implements ProfileListener{ 
	
	private BinningParameters params;
	private PaintableScale scale;
	private int width = 500;
	private int lineWeight=1;
	private boolean lineImageRaster=true;
	private ProfileClusteringHandler clusteringHandler;
	private int numAxisTicks=5;
	private int fontSize=20;
	private Color lineColor =Color.blue;
	private int colorbarHeight=50;
	private boolean addPercentages=false;	
	private Vector<ProfileLinePaintable> linePainters;
	private boolean colorQuantized=false;
	private double [] colorQuantaLimits=null;
	private boolean drawColorBar=true;

	public ProfileLinePanel(BinningParameters bps, PaintableScale s) { 
		params = bps;
		scale = s;
		linePainters = new Vector<ProfileLinePaintable>();
		clusteringHandler = new ProfileClusteringHandler(params);
		
		if(width%params.getNumBins()!=0)
			width = params.getNumBins()*((int)width/params.getNumBins());
		setPreferredSize(new Dimension(width, 300));
	}
	
	public void cluster() { 
		System.out.println("Getting Profiles...");
		Vector<Profile> profs = getAllProfiles();
		System.out.println("Clustering...");
		Vector<Integer> permutation = clusteringHandler.runClustering(profs);
		System.out.println("Reordering...");
		reorder(permutation);
		System.out.println("Done.");
	}
	
	public synchronized Vector<Profile> getAllProfiles() { 
		Vector<Profile> profs = new Vector<Profile>();
		for(int i = 0; i < linePainters.size(); i++) { 
			profs.add(linePainters.get(i).getProfile());
		}
		return profs;
	}

	public synchronized void reorder(Vector<Integer> indices) { 
		Vector<ProfileLinePaintable> newLinePainters = new Vector<ProfileLinePaintable>();
		for(Integer idx : indices) { 
			newLinePainters.add(linePainters.get(idx));
		}
		linePainters = newLinePainters;
		repaint();
	}

	public synchronized void addProfileLinePaintable(ProfileLinePaintable plp) { 
		if(colorQuantized)
			plp.setQuanta(colorQuantaLimits);
		linePainters.add(plp);
		if(linePainters.size() % 10 == 0) { 
			SwingUtilities.invokeLater(new Runnable() { 
				public void run() { 
					updateSize();
				}
			});
			repaint();
		}
	}
	
	public int getPanelWidth(){
		return(width);
	}
	public int getPanelLength(){
		return(colorbarHeight+(linePainters.size()*lineWeight)+lineWeight+1);
	}
	private void updateSize() { 
		setPreferredSize(new Dimension(width, colorbarHeight+(linePainters.size()*lineWeight)+lineWeight+1));
		setSize(new Dimension(width, colorbarHeight+(linePainters.size()*lineWeight)+lineWeight+1));
	}
	
	public void updateFontSize(int size) {
		fontSize = size;
		repaint();
	}
	public void updateLineWeight(int w) {
		lineWeight = w;
		updateSize();
		repaint();
	}
	public void updateColor(Color c) {
		lineColor=c;
		repaint();
	}
	public double getMaxColorVal(){return scale.getMax();}
	public void setMaxColorVal(double v){
		scale.setScale(scale.getMin(), v);
		repaint();
	}
	public double getMinColorVal(){return scale.getMin();}
	public void setMinColorVal(double v){
		scale.setScale(v, scale.getMax());
		repaint();
	}
	public void setDrawColorBar(boolean c){drawColorBar = c;}
	
	public void setLineColorQuanta(double[] q){
		if(q!=null){
			colorQuantaLimits=q;
			colorQuantized=true;
		}
	}
	
	protected void paintComponent(Graphics g) { 
		super.paintComponent(g);
		int w = getWidth(), h = getHeight();
		
		g.setColor(Color.white);
		g.fillRect(0, 0, w, h);

		//Colorbar
		if(drawColorBar)
			drawSiteColorBar((Graphics2D)g, 0, 0);
		
		//Lines
		synchronized(this) { 
			for(int i = 0; i < linePainters.size(); i++) { 
				linePainters.get(i).setColor(lineColor);
				linePainters.get(i).paintItem(g, 0, colorbarHeight+i*lineWeight, w, colorbarHeight+(i*lineWeight)+lineWeight+1);
			}
		}
		
		//Labels, axes, etc 
		g.setColor(Color.DARK_GRAY);
		if(linePainters.size()>0){
			g.drawLine(0, colorbarHeight-1, w, colorbarHeight-1);
			int lastLine = colorbarHeight+(linePainters.size()*lineWeight)+lineWeight+1;
			g.drawLine(0, lastLine, w, lastLine);
		}
		
		//Axis
		g.setFont(new Font("Arial", Font.PLAIN, fontSize));
		FontMetrics metrics = g.getFontMetrics();
		g.setColor(Color.black);
		if(addPercentages){
			if(linePainters.size()>(numAxisTicks*10)){
				for(int i=1; i<=numAxisTicks; i++){
					int l = (100/numAxisTicks)*i;
					String s = String.format("%d%c", l, '%');
					g.drawString(s, w-5-metrics.stringWidth(s), colorbarHeight+(i*((linePainters.size()*lineWeight)/numAxisTicks)));
				}	
			}
		}
	}
	
	public void profileChanged(ProfileEvent p) {
		if(p.getType().equals(ProfileEvent.EventType.ADDED)) { 
			Profile added = p.addedProfile();
			ProfileLinePaintable plp = new ProfileLinePaintable(scale, added);
			scale.updateScale(added.max());
			scale.updateScale(added.min());
			
			//System.out.println(String.format("+Profile: %.3f, %.3f", added.min(), added.max()));
			//System.out.println(String.format("Scale: %.3f, %.3f", scale.getMin(), scale.getMax()));
			
			addProfileLinePaintable(plp);
		}
	}
	public Action createSaveImageAction() { 
	    return new AbstractAction("Save Profile Image...") { 
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
                        saveImage(f, getWidth(), colorbarHeight+(linePainters.size()*lineWeight)+lineWeight+1);
                        //System.out.println("Saved Image [" + sImageWidth + " by " + sImageHeight +  "]");
                    } catch(IOException ie) {
                        ie.printStackTrace(System.err);
                    }
                }
                
            }
        };
	}
	public void saveImage(File f, int w, int h) 
    throws IOException { 
        BufferedImage im = 
            new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        this.paint(g);
        ImageIO.write(im, "png", f);
	}
   
	private void drawSiteColorBar(Graphics2D g2d, int x, int y){
		int cWidth = getWidth();
		int cHeight = colorbarHeight-5;
		
		//Draw colors 
		if(!colorQuantized){
			GradientPaint colorbar = new GradientPaint(x, y, Color.white, x+cWidth, y, lineColor, false);
			g2d.setPaint(colorbar);
			g2d.fillRect(x, y, cWidth, cHeight);
		}else{
			for(int i=0; i<colorQuantaLimits.length; i++){
				int off =(int)x+(i*(int)(cWidth/(double)colorQuantaLimits.length)); 
				g2d.setColor(calcFracColor(lineColor, (double)i/(double)(colorQuantaLimits.length-1)));
				g2d.fillRect(off, y, (int)(cWidth/(double)colorQuantaLimits.length), cHeight);
			}
		}
		
		g2d.setPaint(Color.black);
		g2d.setColor(Color.black);
		
		//Zero tick if needed
		if(scale.getMin()<0){
			int zeroOff = (int)((double)cWidth*(0-scale.getMin())/(scale.getMax()-scale.getMin()));
			g2d.setStroke(new BasicStroke(1.0f));
			g2d.drawLine(zeroOff, 0,zeroOff, cHeight);
		}
		//Draw border
		g2d.setStroke(new BasicStroke(1.0f));
		g2d.drawRect(x, y, cWidth, cHeight);
		
		//Legend
		g2d.setColor(Color.white);
		g2d.setFont(new Font("Ariel", Font.BOLD, 20));
		FontMetrics metrics = g2d.getFontMetrics();
		String maxVal;
		if(scale.getMax()>=1 || scale.getMax()==0)
			maxVal = String.format("%.0f",scale.getMax());
		else
			maxVal = String.format("%.2e",scale.getMax());
		g2d.drawString(maxVal, x+cWidth-(metrics.stringWidth(maxVal)), cHeight-2);
		
	}
	private Color calcFracColor(Color col, double v){
		Color c;
		Color maxColor = col;
		Color minColor = Color.white;
		
		double sVal = v>1 ? 1 : (v<0 ? 0 : v);
		int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
	    int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
	    int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
	    c = new Color(red, green, blue);
		return(c);
	}
}
