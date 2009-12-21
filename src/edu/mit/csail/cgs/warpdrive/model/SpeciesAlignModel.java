package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;

import edu.mit.csail.cgs.datasets.alignments.MultiZAlignRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.MultiZAlignGenerator;
import edu.mit.csail.cgs.warpdrive.components.RegionPanel;

public class SpeciesAlignModel extends WarpModel implements RegionModel, Runnable {

    private static Object theRegionPanels = null;

    private Map<Genome,Map<Genome,MultiZAlignGenerator>> generators;
    private Map<Genome,RegionPanel> regionpanels;
    private Region region;
    private Genome currentGenome;
    private Map<Genome,Region> bestRegion;
    private Map<Genome,ArrayList<MultiZAlignRegion>> alignedRegions;
    private boolean newinput;
    private String alignment;

    public SpeciesAlignModel() {
        generators = new Hashtable<Genome,Map<Genome,MultiZAlignGenerator>>();
        regionpanels = new Hashtable<Genome,RegionPanel>();
        if (theRegionPanels == null) {
            theRegionPanels = regionpanels;
        } else {
            throw new RuntimeException("Sorry, Can't have a second regionpanels");
        }
        bestRegion = new Hashtable<Genome,Region>();
        alignedRegions = new Hashtable<Genome,ArrayList<MultiZAlignRegion>>();
        newinput = false;
        alignment = MultiZAlignGenerator.defaultAlignmentPrefix;
    }

    public void addRegionPanel(RegionPanel p) {
        if (!regionpanels.containsKey(p.getGenome())) {
            p.addModel(this);
            Hashtable<Genome,MultiZAlignGenerator> map = new Hashtable<Genome,MultiZAlignGenerator>();
            for (Genome g : regionpanels.keySet()) {
                map.put(g, new MultiZAlignGenerator(p.getGenome(),
                                                    g));
                Map<Genome,MultiZAlignGenerator> existing = generators.get(g);
                existing.put(p.getGenome(), new MultiZAlignGenerator(g,
                                                                     p.getGenome()));
            }
            generators.put(p.getGenome(),map);
            regionpanels.put(p.getGenome(),p);
        }        
        //        System.err.println("RP Keyset at end of addRegionPanel is " + regionpanels.keySet());
    }

    public void removeRegionPanel(RegionPanel p) {
        p.removeModel(this);
        Genome g = p.getGenome();
        generators.remove(g);
        //        System.err.println("\n\n===============================\nRemoving RP " + p + " which is genome " + g);
        for (Genome o : generators.keySet()) {
            generators.get(o).remove(g);
        }
        regionpanels.remove(g);
        bestRegion.remove(g);
        alignedRegions.remove(g);
    }

    public void setAlignment(String alignment) {
        for (Genome g1 : generators.keySet()) {
            for (Genome g2 : generators.get(g1).keySet()) {
                generators.get(g1).get(g2).setAlignPrefix(alignment);
            }
        }
        this.alignment = alignment;
    }

    public String getAlignment() {return alignment;}

    public void setRegion(Region r) {
        if (newinput) {
            //            System.err.println("Got new input while processing old.  Ignoring");
            return;
        }
        //        System.err.println("SAM going to " + r);
        if (r == bestRegion.get(r.getGenome())) {
            notifyListeners();
        } else {
            region = r;
            newinput = true;
        }
    }

    public synchronized void run() {
        while (keepRunning()) {
            try {
                if (!newinput) {
//                     System.err.println("Waiting... regionpanel keyset is " + regionpanels.keySet());
//                     System.err.println(" regionpanels is " + regionpanels);
//                     if (regionpanels != theRegionPanels) {
//                         throw new RuntimeException("Someone changed regionpanels before wait");
//                     }
                    wait();
//                     System.err.println("Woken... regionpanel keyset is " + regionpanels.keySet());
//                     System.err.println(" regionpanels is " + regionpanels);
//                     if (regionpanels != theRegionPanels) {
//                         throw new RuntimeException("Someone changed regionpanels after wait");
//                     }
                }
            } catch (InterruptedException ex) {}
            if (newinput) {
                try {
                    //                    System.err.println("SAM executing on " + region);
                    Genome g = region.getGenome();
                    int center = (region.getStart() + region.getEnd())/2;
                    int halfsize = (region.getEnd() - region.getStart())/2;
                    currentGenome = g;
                    //                    System.err.println("SAM current Genome is " + g + ".  this is " + this);
                    //                    System.err.println("Current regionpanel keyset is " + regionpanels.keySet());
                    for (Genome o : generators.get(g).keySet()) {
                        if (o.equals(g)) {continue;}
                        //                        System.err.println("  SAM other genome is " + o);
                        RegionPanel p = regionpanels.get(o);
                        if (p == null) {
                            throw new NullPointerException("No RP for Genome " + o );
                        }
                        //                        System.err.println("======  Found RegionPanel for Genome " + o);
                        MultiZAlignRegion best = null;
                        Iterator<MultiZAlignRegion> iter = generators.get(g).get(o).execute(region);
                        ArrayList<MultiZAlignRegion> list = new ArrayList<MultiZAlignRegion>();

                        while (iter.hasNext()) {
                            MultiZAlignRegion mzar = iter.next();
                            list.add(mzar);
                            if (mzar.getStart() <= center &&
                                mzar.getEnd() >= center && 
                                (best == null ||
                                 mzar.getScore() > best.getScore())) {
                                best = mzar;
                            }
                        }            
                        if (best == null) {
                            alignedRegions.put(o,new ArrayList<MultiZAlignRegion>());
                        } else {
                            ArrayList<MultiZAlignRegion> inWindow = new ArrayList<MultiZAlignRegion>();
                            float factor = (float)(center - best.getStart()) / (float)(best.getEnd() - best.getStart());
                            int bestcenter;
                            if (best.getStrand() == '+') {
                                bestcenter = best.getOtherStart() + (int)((best.getOtherEnd() - best.getOtherStart()) * factor);
                            } else {
                                bestcenter = best.getOtherEnd() - (int)((best.getOtherEnd() - best.getOtherStart()) * factor);
                            }

                            Region otherregion = new Region(o,
                                                            best.getOtherChrom(),
                                                            bestcenter - halfsize,
                                                            bestcenter + halfsize);
                            p.setRegion(otherregion);
                            bestRegion.put(o,otherregion);                    
                            for (MultiZAlignRegion r : list) {
                                if (r.getOtherChrom().equals(otherregion.getChrom()) &&
                                    ((r.getOtherStart() >= otherregion.getStart() &&
                                      r.getOtherStart() <= otherregion.getEnd()) || 
                                     (r.getOtherStart() <= otherregion.getStart() &&
                                      r.getOtherEnd() >= otherregion.getStart()))) {
                                    inWindow.add(r);
                                }
                            }
                            alignedRegions.put(o,inWindow);
                        }
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                newinput = false;
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public Region getRegion() {return region;}
    public Genome getCurrentGenome() {return currentGenome;}
    public Set<Genome> getGenomes() {return new HashSet(regionpanels.keySet());}
    public Region getBestRegion(Genome g) {return bestRegion.get(g);}
    public List<MultiZAlignRegion> getAlignedRegions(Genome g) {
        return alignedRegions.get(g);
    }
    /* return the set of alignment version for the current
       set of genomes */
    public Set<String> getAlignments() {
        Set<String> genomes = new HashSet<String>();
        for (Genome g: getGenomes()) {
            genomes.add(g.getVersion());
        }
        //        System.err.println("Looking for alignment versions for " + genomes);
        return MultiZAlignGenerator.getAlignmentVersions(genomes);
    }
}
