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
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Performs hierarchical clustering on a set of motifs and outputs the clusters
 * along with PNGs of the clusters.
 */

public class ClusterMotifs {    
    private WMComparator comparator;
    private ClusteringMethod<WeightMatrix> method;
    private ClusterRepresentative<WeightMatrix> rep;
    private int minClusterSize;

    public static void main(String args[]) {
        try {
            boolean normalize = Args.parseFlags(args).contains("normalize");
            int compareLength = Args.parseInteger(args,"compareLength",-1);
            double maxDistance = Args.parseDouble(args,"maxDistance",3.0);
            String pictureDirectory = Args.parseString(args,"pictures",null);
            int minClusterSize = Args.parseInteger(args,"minclustersize",2);
            Collection<WeightMatrix> allmatrices = Args.parseWeightMatrices(args);

            
            MarkovBackgroundModel bgModel = null;
            String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
            BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                                  1,
                                                                                  "MARKOV",
                                                                                  Args.parseGenome(args).cdr().getDBID());
            if (md != null) {
                bgModel = BackgroundModelLoader.getMarkovModel(md);
            } else {
                System.err.println("Couldn't get metadata for " + bgmodelname);
            }


            for (WeightMatrix m : allmatrices) {
                m.toFrequency(bgModel);
            }

            ClusterMotifs cm = new ClusterMotifs(maxDistance,compareLength,normalize, minClusterSize);
            Collection<Cluster<WeightMatrix>> output = cm.cluster(allmatrices);
            int count = 0;
            Cluster sorted[] = new Cluster[output.size()];            
            for (Cluster<WeightMatrix> cluster : output) {
                sorted[count++] = cluster;
            }
            Arrays.sort(sorted,new ClusterSizeComparator());
            for (count = 0; count < sorted.length; count++) {
                if ((sorted[count].size() >= minClusterSize)) {
                    System.out.println(count + "=======\n" + cm.printCluster(sorted[count]) + "\n");
                    if (pictureDirectory != null) {
                        cm.drawCluster(sorted[count].getElements(),pictureDirectory +"/cluster" +count + ".png");
                    }
                }
            }
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
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

    public ClusterMotifs(double maxdist, int compareLength, boolean normalize, int minClusterSize) {
        comparator = new WMDistanceComparator(normalize,compareLength);
        rep = new WMMinAvgDistanceRep(comparator);
        method = new HierarchicalClustering<WeightMatrix>(rep,comparator);
        ((HierarchicalClustering)method).setMaxDistanceToAccept(maxdist);
        this.minClusterSize = minClusterSize;
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
        WeightMatrix elements[] = new WeightMatrix[wms.size()];
        int i = 0; 
        for (WeightMatrix m : wms) {
            if (m == null) {
                throw new NullPointerException("null weight matrix in cluster " + cluster);
            }
            if (m.name == null) {
                throw new NullPointerException("null name at " + (i-1) + " in cluster " + cluster);
            }
            elements[i++] = m;
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
    public static void drawCluster(Collection<WeightMatrix> cluster, String fname) {
        final int pixwidth = 800;
        final int pixheight = 200;
        final int maxcols = 4;
        drawCluster(cluster,fname,pixwidth,pixheight,maxcols);
    }
    public static void drawCluster(Collection<WeightMatrix> cluster, String fname, int pixwidth, int pixheight, int maxcols) {
        File f = new File(fname);
        int rows,cols;
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
        Iterator<WeightMatrix> matrices = cluster.iterator();
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


