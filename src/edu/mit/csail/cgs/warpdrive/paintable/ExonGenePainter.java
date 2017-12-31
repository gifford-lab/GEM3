package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.sql.SQLException;
import java.util.*;

import javax.swing.JToolTip;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.ExonicGene;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.components.GOAnnotationPanel;
import edu.mit.csail.cgs.warpdrive.model.*;

public class ExonGenePainter extends RegionPaintable {
    
    private GeneModel model;
    private NonOverlappingLayout<Gene> layout;

    private FunctionLoader funcLoader=null;
    private GOAnnotationPanel.Frame goFrame;    
    private GeneProperties props;
    private DynamicAttribute attrib;
    private double htRat, wdRat;
    private Vector<Gene> genes;    
    
    public ExonGenePainter(GeneModel model) {
        super();
        layout = new NonOverlappingLayout<Gene>();
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        htRat = .03;
        wdRat = .03;
        model.addEventListener(this);
        goFrame = null;
        props = new GeneProperties();
        initLabels();
/* I don't think anyone uses the GO annotations now. 
        try {
            funcLoader = new GOFunctionLoader("go_200904");			
        } catch (SQLException e) {
            e.printStackTrace();
            if(funcLoader != null) { funcLoader.close(); }
            funcLoader = null;
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            if(funcLoader != null) { funcLoader.close(); }
            funcLoader = null;
        }
*/
    }

    public GeneProperties getProperties() {
        return props;
    }

    public void clickedOnItem(ActionEvent e) {
        String geneName = e.getActionCommand();
        
        if(funcLoader != null) { 
            if(goFrame == null) { 
                try {
                    GOAnnotationPanel panel = new GOAnnotationPanel(funcLoader, getRegion().getGenome());
                    panel.setID(geneName);
                    panel.setVersion(getRegion().getGenome().getSpecies());
                    panel.annotate();
                    
                    goFrame = new GOAnnotationPanel.Frame(panel);
                } catch (SQLException e1) {
                    e1.printStackTrace();
                }
                
            } else { 
                GOAnnotationPanel panel = goFrame.getPanel();
                panel.setID(geneName);
                panel.annotate();
                goFrame.setVisible(true);
            }
        }
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
        if(funcLoader != null)
        	funcLoader.close();
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
            genes = null;
            setLayoutGenes();
            notifyListeners();
        }
    }    

    public int getMaxVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }

    private void setLayoutGenes() {
        if (canPaint() && genes == null) {
            Iterator<Gene> itr = model.getResults();
            genes = new Vector<Gene>();
            while(itr.hasNext()) { genes.add(itr.next()); }
            layout.setRegions(genes); 
        }
    }
    
    public int getNumTracks() { return layout.getNumTracks(); }
    
    private static int xcoord(int base, int x1, int baseStart, double scale) { 
        return x1 + (int)Math.round((double)(base - baseStart) * scale);
    }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        if (!canPaint()) {
            return;
        }
        boolean drawgenenames = props.DrawGeneNames;
        boolean drawallgenenames = props.DrawAllGeneNames;
               
        int w = x2 - x1, h = y2 - y1;
 
        int numTracks = Math.max(1, layout.getNumTracks());
        int trackHeight = Math.max(2, h/(numTracks*2));
        int halfTrackHeight = Math.max(1, trackHeight/2);
        
        Region region = model.getRegion();
        int rs = region.getStart(), re = region.getEnd();
        int rw = re - rs + 1;
        double xScale = (double)w / (double)rw;
        int[] a = new int[7];
        int[] b = new int[7];
        
        clearLabels();
        boolean drewAnything = false;
        
        Font oldFont = g.getFont();
        g.setFont(attrib.getRegionLabelFont(w,h));
        FontMetrics fontmetrics = g.getFontMetrics();
        //--------------------------------------------------
        
        HashSet<String> labels = new HashSet<String>();

        for(Gene gene : genes) {
            int track = 0;
            
            if(!layout.hasTrack(gene)) { 
                System.err.println("No track assigned to gene: " + gene.getName());
            } else { 
                track = layout.getTrack(gene);
            }

            int gy1 = y1 + (2 * trackHeight * track);
            int gy2 = gy1 + trackHeight;
            int texty = gy2 + trackHeight;
            int gmy = gy1 + halfTrackHeight;
            
            int geneHeight = Math.max(2, (int)Math.floor((double)trackHeight * 0.80));
            int halfHeight = trackHeight / 2;
            int halfGeneHeight = geneHeight / 2;
            
            int gtop = gmy - (halfGeneHeight/2), gbottom = gtop + (geneHeight/2);
            int rectheight = gbottom - gtop;
            
            int geneStart = gene.getStart(), geneEnd = gene.getEnd();
            boolean strand = gene.getStrand() == '+';
            
            int gx1 = xcoord(geneStart, x1, rs, xScale);
            int gx2 = xcoord(geneEnd, x1, rs, xScale);
            int gleft = Math.max(x1, gx1), gright = Math.min(x2, gx2);

            g.setColor(Color.black);
            g.drawLine(gx1, gmy, gx2, gmy);
            arrangeArrow(a, b, strand, trackHeight, gx1, gx2, gmy);
            g.drawPolyline(a, b, 7);
            
            if(gene instanceof ExonicGene) { 
                ExonicGene exonGene = (ExonicGene)gene;
                
                Iterator<Region> exons = exonGene.getExons();
                while(exons.hasNext()) { 
                    Region exon = exons.next();
                    int ex1 = xcoord(exon.getStart(), x1, rs, xScale);
                    int ex2 = xcoord(exon.getEnd(), x1, rs, xScale);
                    int eleft = Math.max(x1, ex1);
                    int eright = Math.min(x2, ex2);

                    int rectwidth = eright - eleft + 1;

                    g.setColor(Color.pink);
                    g.fillRect(eleft, gtop, rectwidth, gbottom - gtop);
                    g.setColor(Color.black);
                    g.drawRect(eleft, gtop, rectwidth, gbottom - gtop);

                    drewAnything = true;
                }
                
            } else {
                int rectwidth = gright - gleft + 1;

                g.setColor(Color.pink);
                g.fillRect(gleft, gtop, rectwidth, gbottom - gtop);
                g.setColor(Color.black);
                g.drawRect(gleft, gtop, rectwidth, gbottom - gtop);

                drewAnything = true;
            }
            
                
            boolean alwaysDrawNames = props.AlwaysDrawNames;
            // Somewhat prettier gene-name output

            int fontsize = g.getFont().getSize();
            int nx = Math.max(x1 + 3, gx1 + 3);  // gotta do the Math.max(), to make sure the name is on the screen.

            //int ny = gmy + (halfGeneHeight / 2) - 2;
            //int ny = texty;
            int ny = gmy + (halfGeneHeight/2) + Math.min(fontsize, trackHeight) + 1;
            
            int rectwidth = gright - gleft + 1;

            ArrayList<String> aliases = new ArrayList<String>();
            labels.clear();            

            String first;
            if (gene.getName().endsWith("Rik")) {
                first = gene.getID();
                aliases.add(gene.getName());
            } else {
                first = gene.getName();
                aliases.add(gene.getID());
            }
            aliases.addAll(gene.getAliases());
            String todraw = null;
            addLabel(gleft,gy1,rectwidth,trackHeight*2,first);

            labels.add(first);
            if(alwaysDrawNames || first.length() * fontsize < rectwidth) {
                todraw = first;
            }

            g.setColor(Color.black);
            for (String s : aliases) {
            	if(!s.endsWith("Rik")) { 
            		String newtodraw = todraw + ", " + s;
            		if (drawallgenenames && todraw != null && 
            				fontmetrics.charsWidth((newtodraw).toCharArray(),0,newtodraw.length()) < rectwidth) {
            			todraw = newtodraw;
            		}

            		if(!labels.contains(s)) { 
            			addLabel(gleft,gy1,rectwidth,trackHeight*2,s);
            			labels.add(s);
            		}
            	}
            }

            if (todraw != null && drawgenenames) {
                Rectangle2D textrect = fontmetrics.getStringBounds(todraw, g);
                int diff = (gright-gleft)/2 - (int)Math.round(textrect.getWidth()) / 2;
                g.drawString(todraw, nx + diff, ny);
            }
        }
        
        g.setFont(oldFont);
    }
    
    private void arrangeArrow(int[] a, int[] b, boolean strand, int geneHeight, int gx1, int gx2, int my) { 
        int gxw = gx2 - gx1;
        double arrowHt = htRat * geneHeight;
        double arrowWd = wdRat * gxw;
        int a1, a2, a3;
        
        if(arrowWd > arrowHt) { 
            arrowWd = arrowHt; 
        }
        
        // forward arrow
        if(strand) {
            int startX = gx1;
            a1 = startX; 
            a2 = (int) Math.round(startX + (arrowWd * 8)); 
            a3 = (int) Math.round(startX + (arrowWd * 12)); 
            
        } else { 
            // backward arrow
            int startX = gx2+1;
            a1 = startX; 
            a2 = (int) Math.round(startX - (arrowWd * 8)); 
            a3 = (int) Math.round(startX - (arrowWd * 12)); 
        }
        
        a[0] = a1;
        a[1] = a1;
        a[2] = a2;
        a[3] = a2;
        a[4] = a3;
        a[5] = a2;
        a[6] = a2;

        int b1 = (int) Math.round(my);
        int b2 = (int) Math.round(my - (arrowHt * 13));
        int b3 = (int) Math.round(my - (arrowHt * 10));
        int b4 = (int) Math.round(my - (arrowHt * 16));
        
        b[0] = b1;
        b[1] = b2;
        b[2] = b2;
        b[3] = b3;
        b[4] = b2;
        b[5] = b4;
        b[6] = b2;
    }

}
