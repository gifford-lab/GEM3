/*
 * Created on Oct 19, 2007
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.motifs.ReverseComplement;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.viz.utils.FileChooser;

public class PREScanner implements Mapper<Region,Double[]> {

    public static void main(String[] args) { 
    	FileChooser chooser = new FileChooser(null);
        int width = args.length > 0 ? Integer.parseInt(args[0]) : 220;
        File input = args.length > 1 ? new File(args[1]) : chooser.choose();
        
        try {
            PRELoader loader = new PRELoader(width, input);
            PREScanner scanner = loader.getScanner();
            
            String line;
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            Pattern rp = Pattern.compile("(.*):(\\d+)-(\\d+)");
            Genome g = Organism.findGenome("mm8");
            
            System.out.print(">"); System.out.flush();
            while((line = br.readLine()) != null) { 
                if((line = line.trim()).length() > 0) { 
                    Matcher m = rp.matcher(line);
                    if(m.matches()) { 
                        String chrom = m.group(1);
                        int start = Integer.parseInt(m.group(2));
                        int end = Integer.parseInt(m.group(3));
                        Region r = new Region(g, chrom, start, end);
                        
                        Double[] scores = scanner.scoreRegion(r);
                        double max = 0.0;
                        for(int i = 0; i < scores.length; i++) { 
                            System.out.print(String.format("%.1f ", scores[i]));
                            max = Math.max(max, scores[i]);
                        }
                        System.out.println(String.format("\nMax Score: %.1f", max));
                    }
                }
                System.out.print(">"); System.out.flush();
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
    
    private int width;
    private Map<String,MotifScanner> scanners;
    private Map<String,Double> singleWeights;
    private Map<String,Map<String,Double>> pairWeights;
    
    public PREScanner(int w, Vector<MotifScanner> scs, Map<String,Double> sws, Map<String,Map<String,Double>> pws) { 
        width = w;
        scanners = new HashMap<String,MotifScanner>();
        
        for(MotifScanner ms : scs) { 
            scanners.put(ms.getName(), ms);
        }
        
        singleWeights = new HashMap<String,Double>(sws);
        pairWeights = new HashMap<String,Map<String,Double>>(pws);
    }
    
    public PREScanner(int width, MotifScanner[] scs, double[] sws, double[][] pws) { 
        this.width = width;
        
        scanners = new HashMap<String,MotifScanner>();
        singleWeights = new HashMap<String,Double>();
        pairWeights = new HashMap<String,Map<String,Double>>();
        
        for(int i = 0; i < scs.length; i++) { 
            scanners.put(scs[i].getName(), scs[i]);
            singleWeights.put(scs[i].getName(), sws[i]);
            pairWeights.put(scs[i].getName(), new HashMap<String,Double>());

            for(int j = 0; j <= i; j++) { 
                pairWeights.get(scs[j].getName()).put(scs[i].getName(), pws[i][j]);
                pairWeights.get(scs[i].getName()).put(scs[j].getName(), pws[i][j]);
            }
        }
    }
    
    public double scoreRegionAtPoint(Point p) { 
    	int start = p.getLocation()-width/2;
    	int end = p.getLocation()+width/2;
    	Region window = new Region(p.getGenome(), p.getChrom(), start, end);

        RunningMotifState ms = 
            new RunningMotifState(scanners.keySet(), singleWeights, pairWeights);
        
        for(String n : scanners.keySet()) {
            Site[] array = scanners.get(n).findSites(window);
            for(int i = 0; i < array.length; i++) { 
            	Site s = array[i];
            	ChangePoint cp = new ChangePoint(n, s.getStart(), 1);
            	ms.updateState(cp);
            }
        }

        return ms.calculateScore();
    }
    
    public Double[] scoreRegion(Region r) {
        int w2 = width/2;
        Double[] scores = new Double[r.getEnd()-r.getStart()+1];
        Region wr = new Region(r.getGenome(), r.getChrom(), 
                r.getStart()-w2, r.getEnd()+w2);

        Vector<ChangePoint> cpvec = new Vector<ChangePoint>();
        Map<String,Site[]> sites = new HashMap<String,Site[]>();
        
        for(String n : scanners.keySet()) {
            Site[] array = scanners.get(n).findSites(wr);
            sites.put(n, array);
            for(int i = 0; i < array.length; i++) { 
                cpvec.add(new ChangePoint(n, array[i].getStart() - w2, 1));
                cpvec.add(new ChangePoint(n, array[i].getEnd() + w2, -1));
            }
        }
        
        ChangePoint[] cparray = cpvec.toArray(new ChangePoint[cpvec.size()]);
        Arrays.sort(cparray);
        
        RunningMotifState ms = 
            new RunningMotifState(scanners.keySet(), singleWeights, pairWeights);
        
        int cpi = 0;
        for(int i = 0; i < scores.length; i++) {
            int offset = r.getStart()+i;
            while(cpi < cparray.length && cparray[cpi].getChangeLocation() <= offset) { 
                ms.updateState(cparray[cpi++]);
            }
            scores[i] = ms.calculateScore();
        }
        
        return scores;
    }

    public Double[] execute(Region a) {
        return scoreRegion(a);
    }
}

class RunningMotifState {
    
    private String[] motifArray;
    private Map<String,int[]> runningCounts;
    private Map<String,Double> singleWeights; 
    private Map<String,Map<String,Double>> pairWeights;
    private double currentWeight;
    
    public RunningMotifState(Set<String> motifNames, 
            Map<String,Double> sw, Map<String,Map<String,Double>> pw) { 
        runningCounts = new HashMap<String,int[]>();
        for(String n : motifNames) { 
            runningCounts.put(n, new int[1]);
            runningCounts.get(n)[0] = 0;
        }
        
        motifArray = motifNames.toArray(new String[motifNames.size()]);
        Arrays.sort(motifArray);
        
        singleWeights = new HashMap<String,Double>(sw);
        pairWeights = new HashMap<String,Map<String,Double>>(pw);
        
        currentWeight = 0.0;
    }
    
    public double getScore() { return currentWeight; }
    
    public double calculateScore() {
        double weight = 0.0;
        
        for(int i = 0; i < motifArray.length; i++) {
            String mi = motifArray[i];
            
            weight += getMotifCount(mi) * singleWeights.get(mi);
            
            /*
            if(hasCurrentMotif(mi)) { 
                weight += singleWeights.get(mi);
            }
            */
            
            for(int j = 0; j <= i; j++) {
                String mj = motifArray[j];
                if(hasCurrentMotifPair(mi, mj)) { 
                    weight += pairWeights.get(mi).get(mj);
                } else { 
                	weight -= pairWeights.get(mi).get(mj);
                }
            }
        }
        
        return weight;
    }

    public int getCount(String m) { return runningCounts.get(m)[0]; }
    public Set<String> getPossibleMotifs() { return runningCounts.keySet(); }
    
    public boolean hasCurrentMotif(String m) { return runningCounts.get(m)[0] > 0; }
    public int getMotifCount(String m) { return runningCounts.get(m)[0]; }
    
    public boolean hasCurrentMotifPair(String m1, String m2) { 
        return Math.min(1, getCount(m1)) + Math.min(1, getCount(m2)) >= 2; 
    }
    
    public Set<String> getCurrentMotifs() { 
        HashSet<String> ms = new HashSet<String>();
        for(String m : getPossibleMotifs()) { 
            if(getCount(m) > 0) { 
                ms.add(m);
            }
        }
        return ms;
    }
    
    public void updateState(ChangePoint cp) {
        int currentCount = runningCounts.get(cp.getName())[0];
        if(currentCount + cp.getCountDiff() < 0) { throw new IllegalArgumentException(); }
        
        boolean added = currentCount == 0 && cp.getCountDiff() > 0;
        boolean removed = currentCount > 0 && currentCount + cp.getCountDiff() == 0;
        
        if(removed) { 
            //currentWeight -= singleWeights.get(cp.getName());
            for(String m : getPossibleMotifs()) { 
                if(hasCurrentMotifPair(cp.getName(), m)) { 
                    currentWeight -= pairWeights.get(cp.getName()).get(m);
                }
            }            
        }
        
        currentWeight += (cp.getCountDiff() * singleWeights.get(cp.getName()));        
        runningCounts.get(cp.getName())[0] += cp.getCountDiff();
        
        if(added) { 
            //currentWeight += singleWeights.get(cp.getName());
            for(String m : getPossibleMotifs()) { 
                if(hasCurrentMotifPair(cp.getName(), m)) { 
                    currentWeight += pairWeights.get(cp.getName()).get(m);
                }
            }
        }
    }
}

class ChangePoint implements Comparable<ChangePoint> { 
    
    public ChangePoint(String n, int cloc, int cd) { 
        name = n;
        changeLocation = cloc;
        countDiff = cd;
    }
    
    private int changeLocation;
    private String name;
    private int countDiff;
    
    public int getChangeLocation() { return changeLocation; }
    public int getCountDiff() { return countDiff; }
    public String getName() { return name; }
    
    public boolean equals(Object o) {
        if(!(o instanceof ChangePoint)) { return false; }
        ChangePoint cp = (ChangePoint)o;
        if(!name.equals(cp.name)) { return false;  }
        if(changeLocation != cp.changeLocation) { return false; }
        return true;
    }
    
    public int compareTo(ChangePoint cp) { 
        if(changeLocation < cp.changeLocation) { return -1; }
        if(changeLocation > cp.changeLocation) { return 1; }
        return 0;
    }
}

class Site extends StrandedRegion {
    
    private String motifName;
    private String sequence;
    
    public Site(Genome g, String c, int start, int end, char str, 
            String motifName, String seq) {
        super(g, c, start, end, str);
        this.motifName = motifName;
        sequence = seq;
    }
    
    public String getMotifName() { return motifName; }
    public String getSequence() { return sequence; }
}

interface MotifScanner { 
    public int getWidth();
    public String getName();
    public Site[] findSites(Region r);
}

class RegexMotifScanner implements MotifScanner {
    
    public static SequenceGenerator seqgen;
    public static ReverseComplement revcomp;
    
    static {
        seqgen = new SequenceGenerator();
        revcomp = new ReverseComplement();
    }
    
    private String name;
    private Pattern pattern;
    private int width;
    
    public RegexMotifScanner(String n, String p) { 
        name = n;
        pattern = Pattern.compile(p);
        width = p.length();
    }

    public Site[] findSites(Region r) {
        Vector<Site> sites = new Vector<Site>();
        
        String fstrand = seqgen.execute(r).toUpperCase();
        sites.addAll(findStringSites(r, fstrand, '+'));
        
        String rstrand = revcomp.execute(fstrand);
        sites.addAll(findStringSites(r, rstrand, '-'));
        
        return sites.toArray(new Site[sites.size()]);
    }
    
    public Collection<Site> findStringSites(Region r, String str, char strand) { 
        LinkedList<Site> sites = new LinkedList<Site>();
        
        Matcher m = pattern.matcher(str);
        while(m.find()) { 
            int offsetstart = m.start(), offsetend = m.end();
            int start = r.getStart() + offsetstart;
            int end = r.getStart() + offsetend;
            Site s = new Site(r.getGenome(), r.getChrom(), start, end, strand, 
                    name, m.group(0));
            sites.addLast(s);
        }
        
        return sites;
    }

    public String getName() {
        return name;
    }

    public int getWidth() {
        return width;
    } 
}
