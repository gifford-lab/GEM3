package edu.mit.csail.cgs.datasets.motifs;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.font.*;


public class WeightMatrixPainter {
	public final static int X_MARGIN = 5;
	public final static int Y_MARGIN = 2;
	public final static int YLABEL_SIZE = 12;
    /* paints a representation of a weight matrix wherein the height of the letters roughly indicates
       their relative probability at each base.  The image is painted in g in the
       bounding box defined by the upper left corner x1,y1 and lower right corner x2,y2 */
    public void paint(WeightMatrix wm, Graphics g, int x1, int y1, int x2, int y2) {

        wm.toLogOdds();

        int w = x2 - x1;
        int h = y2 - y1;
        Font labelFont = new Font("Arial",Font.PLAIN,YLABEL_SIZE);
        g.setFont(labelFont);
        FontMetrics fontmetrics = g.getFontMetrics();
        String label = wm.toString();
        LineMetrics linemetrics = fontmetrics.getLineMetrics(label,g);
        g.setColor(Color.BLACK);
        g.drawString(label,x1+X_MARGIN + w/2 - fontmetrics.charsWidth(label.toCharArray(),0,label.length()) / 2,y2-Y_MARGIN);
        int labelHeight = fontmetrics.getHeight() + Y_MARGIN;

        Font baseFont = new Font("Arial",Font.BOLD, h );
        int pixelsPerLetter = (w-X_MARGIN*2) / wm.length();
        double xfactor = ((float)pixelsPerLetter) / baseFont.getSize() * 1.3;
        //        System.err.println("Xfactor is " + xfactor + " and base sizeis " + baseFont.getSize());
        double vals[] = new double[4];
        for (int pos = 0; pos < wm.length(); pos++) {
            Character[] letters = WeightMatrix.getLetterOrder(wm,pos);
            double total = 0;
            int ypos = y2 - labelHeight;
            for (int j = 3; j >= 0; j--) {
                vals[j] = Math.exp(wm.matrix[pos][letters[j]]);
                total += vals[j];
            }
            
            double bits = 0;
            for (int j = 3; j >= 0; j--) {
                vals[j] = vals[j] / total;
                bits -= vals[j] * (Math.log(vals[j]) / Math.log(2));
            }
            double totalscale = (2.0 - bits) / 2.0;
            for (int j = 3; j >= 0; j--) {
                double val = vals[j];
                AffineTransform transform = new AffineTransform(xfactor,0,0,val * totalscale,0,0);
                Font thisFont = baseFont.deriveFont(transform);
                g.setFont(thisFont);
                int offset = 0;
                if (letters[j] == 'A') {
                    g.setColor(Color.GREEN);
                } else if (letters[j] == 'C') {
                    g.setColor(Color.BLUE);
//                    offset = (int)(thisFont.getSize()*0.04);
                } else if (letters[j] == 'G') {
                    g.setColor(Color.ORANGE);
                    offset = -(int)(pixelsPerLetter*0.05);
                }  else if (letters[j] == 'T') {
                    g.setColor(Color.RED);
                    offset = (int)(pixelsPerLetter*0.05);
                }
                g.drawString(letters[j].toString(),x1+X_MARGIN +offset+ pos * pixelsPerLetter,ypos);
                ypos -= thisFont.getSize() * val * totalscale;
            }                
        }
    }

}
