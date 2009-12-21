package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.sql.*;
import java.io.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.clustering.*;
import edu.mit.csail.cgs.clustering.hierarchical.*;

public class ClusterMotifs {

    

    
    private WMComparator comparator;
    private ClusteringMethod<WeightMatrix> method;
    private ClusterRepresentative<WeightMatrix> rep;

    public static void main(String args[]) {
        try {
            boolean normalize = false;
            int compareLength = -1;
            double maxDistance = 3;
            String pictureDirectory = null;
            for (int i = 0; i < args.length;i++) {
                if (args[i].equals("--normalize")) {
                    normalize = true;
                }
                if (args[i].equals("--maxDistance")) {
                    maxDistance = Double.parseDouble(args[++i]);
                }
                if (args[i].equals("--compareLength")) {
                    compareLength = Integer.parseInt(args[++i]);
                }
                if (args[i].equals("--pictures")) {
                    pictureDirectory = args[++i];
                }
            }

            String getmatrices = "select id from weightmatrix";   
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement stmt = cxn.prepareStatement(getmatrices);            
            ResultSet rs = stmt.executeQuery();
            ArrayList<WeightMatrix> allmatrices = new ArrayList<WeightMatrix>();
            while (rs.next()) {
                int wmid = rs.getInt(1);
                WeightMatrix m = WeightMatrix.getWeightMatrix(wmid);
                m.toFrequency();
                if (m.length() > 5) {
                    allmatrices.add(m);
                }
            }
            rs.close();
            stmt.close();

            ClusterMotifs cm = new ClusterMotifs(maxDistance,compareLength,normalize);
            Collection<Cluster<WeightMatrix>> output = cm.cluster(allmatrices);
            int count = 0;
            Cluster sorted[] = new Cluster[output.size()];            
            for (Cluster<WeightMatrix> cluster : output) {
                sorted[count++] = cluster;
            }
            Arrays.sort(sorted,new ClusterSizeComparator());
            for (count = 0; count < sorted.length; count++) {
                if ((sorted[count].size() > 1)) {
                    System.out.println(count + "=======\n" + cm.printCluster(sorted[count]) + "\n");
                    if (pictureDirectory != null) {
                        cm.drawCluster(sorted[count],pictureDirectory +"/cluster" +count + ".png");
                    }
                }
            }
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        }
    }

    public static boolean hasType(Cluster<WeightMatrix> cluster, String type) {
        for (WeightMatrix wm : cluster.getElements()) {
            if (wm.type.equals(type)) {
                return true;
            }
        }
        return false;
    }

    public ClusterMotifs(double maxdist, int compareLength, boolean normalize) {
        comparator = new WMDistanceComparator(normalize,compareLength);
        rep = new WMMinAvgDistanceRep(comparator);
        method = new HierarchicalClustering<WeightMatrix>(rep,comparator);
        ((HierarchicalClustering)method).setMaxDistanceToAccept(maxdist);
    }

    public ClusterMotifs(WMComparator comp,
                         ClusterRepresentative<WeightMatrix> rep,
                         ClusteringMethod<WeightMatrix> method) {
        this.comparator = comp;
        this.rep = rep;
        this.method = method;
    }

    public Collection<Cluster<WeightMatrix>> cluster(Collection<WeightMatrix> elements) {
        return method.clusterElements(elements);
    }
    
    public String printCluster(Cluster<WeightMatrix> cluster) {
        StringBuffer b = new StringBuffer();
        Set<WeightMatrix> wms = cluster.getElements();
        Iterator<WeightMatrix> iter = wms.iterator();
        WeightMatrix elements[] = new WeightMatrix[cluster.size()];
        int i = 0; 
        while (iter.hasNext()) {
            elements[i++] = iter.next();
        }
        Arrays.sort(elements, new WeightMatrixSorter());
        b.append("Matrices : ");
        for (i = 0; i < elements.length; i++) {            
            if (i > 0) {b.append(", ");}
            b.append(elements[i].name + "(" + elements[i].version + ")");
        }
        b.append("\n");
        b.append(WeightMatrix.printMatrixLetters(rep.getRepresentative(cluster)));
        return b.toString();
    } 
    public void drawCluster(Cluster<WeightMatrix> cluster, String fname) {
        File f = new File(fname);
        final int pixwidth = 800;
        final int pixheight = 200;
        int rows,cols;
        final int maxcols = 4;
        double squareness[] = new double[maxcols];
        for (cols = 1; cols <= maxcols; cols++) {
            rows = (int)Math.ceil(cluster.size() / (float)cols);        
            double sq = Math.max((cols * pixwidth) / (rows * pixheight), (rows*pixheight) / (cols * pixwidth));
            //            System.err.println("size " + cluster.size() + " cols=" + cols + "  => " + sq);
            squareness[cols - 1] = sq;
        }
        double min = 100000; int minind = -1;
        for (cols = 0; cols < maxcols; cols++) {
            if (squareness[cols] < min) {
                min = squareness[cols];
                minind = cols;
            }
        }
        cols = minind + 1;
        rows = (int)Math.ceil(cluster.size() / (float)cols);        
        //        System.err.println("Trying to use cols " + cols + ", rows " + rows + " for size " + cluster.size());

        BufferedImage im = new BufferedImage(pixwidth * cols, pixheight * rows,BufferedImage.TYPE_INT_RGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        WeightMatrixPainter wmp = new WeightMatrixPainter();
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,pixwidth * cols, pixheight * rows);
        Iterator<WeightMatrix> matrices = cluster.getElements().iterator();
        for (int i = 0; i < cluster.size(); i++) {
            int c = i % cols;
            int r = (i / cols);
            //            System.err.println(" i -> " + c + "," + r);
            wmp.paint(matrices.next(),g2,c * pixwidth,r * pixheight,(c + 1 ) * pixwidth,(r + 1) * pixheight);
        }
        try {
            ImageIO.write(im,"png",f);
        }  catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
class WeightMatrixSorter implements Comparator<WeightMatrix> {
    public int compare(WeightMatrix a, WeightMatrix b) {
        return a.name.compareTo(b.name);
    }
}

class ClusterSizeComparator implements Comparator<Cluster> {
    public int compare(Cluster a, Cluster b) {
        return a.size() - b.size();
    }
}


