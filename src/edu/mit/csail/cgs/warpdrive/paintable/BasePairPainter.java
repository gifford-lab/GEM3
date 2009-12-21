package edu.mit.csail.cgs.warpdrive.paintable;
import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionMapperModel;

public class BasePairPainter extends RegionPaintable {

    private RegionMapperModel<String> model;
    private Color A = Color.GREEN, C = Color.BLUE, G = Color.ORANGE, T = Color.RED;
    private BasePairPainterProperties props;

    public BasePairPainter (RegionMapperModel<String> model) {
        super();
        this.model = model;
        model.addEventListener(this);
        props = new BasePairPainterProperties();
    }
    
    public BasePairPainterProperties getProperties() {
        return props;
    }

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public synchronized void eventRegistered(EventObject e) {
        if (e.getSource() == model &&
            model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }    

    public int getMaxVertSpace() { 
        return props.FontSize * 2;
    }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        if (!canPaint()) {
            return;
        }
        int w = x2 - x1;
        double h = y2 - y1;
        Region region = getRegion();
        int regionstart = region.getStart(), regionend = region.getEnd();
        int regionwidth = regionend - regionstart;
        int fontsize = props.FontSize;
        String wholestring = model.getResults();
        if (wholestring == null) {
            return;
        }
        char[] chars = wholestring.toUpperCase().toCharArray();
        Font oldFont = g.getFont();
        Font newFont = new Font("Arial",Font.PLAIN,fontsize);
        g.setFont(newFont);
        if (chars.length > w) {
            // need to draw some summary that combines multiple bases into single pixels

        } else if (chars.length * fontsize < w) {
            // can draw one character per 
            for (int i = 0; i < chars.length; i++) {
                if (chars[i] == 'A') {
                    g.setColor(A);
                } else if (chars[i] == 'C') {
                    g.setColor(C);
                } else if (chars[i] == 'T') {
                    g.setColor(T);
                } else if (chars[i] == 'G') {
                    g.setColor(G);
                } else {
                    g.setColor(Color.GRAY);
                }
                int x = getXPos(regionstart + i,
                                regionstart,
                                regionend,
                                x1,
                                x2);
                g.drawString(Character.toString(chars[i]),x, y2);
            }
        } else if (chars.length * fontsize < 2*w) {
            // can draw one character per 
            for (int i = 0; i < chars.length; i++) {
                if (chars[i] == 'A') {
                    g.setColor(A);
                } else if (chars[i] == 'C') {
                    g.setColor(C);
                } else if (chars[i] == 'T') {
                    g.setColor(T);
                } else if (chars[i] == 'G') {
                    g.setColor(G);
                } else {
                    g.setColor(Color.GRAY);
                }
                int x = getXPos(regionstart + i,
                                regionstart,
                                regionend,
                                x1,
                                x2) + (i % 2);
                g.drawString(Character.toString(chars[i]),x, y2 - (i % 2) * fontsize);
            }
        } else {
            // draw color block representation
            for (int i = 0; i < chars.length; i++) {
                if (chars[i] == 'A') {
                    g.setColor(A);
                } else if (chars[i] == 'C') {
                    g.setColor(C);
                } else if (chars[i] == 'T') {
                    g.setColor(T);
                } else if (chars[i] == 'G') {
                    g.setColor(G);
                } else {
                    g.setColor(Color.GRAY);
                }
                int xleft = getXPos(regionstart + i,
                                    regionstart,
                                    regionend,
                                    x1,
                                    x2);
                int xright = getXPos(regionstart + i + 1,
                                     regionstart,
                                     regionend,
                                     x1,
                                     x2);
                g.fillRect(xleft,y1,xright - xleft, y2-y1);
            }
        }
        g.setFont(oldFont);
    }
}
