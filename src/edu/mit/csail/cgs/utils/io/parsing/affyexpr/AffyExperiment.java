/*
 * Created on Sep 6, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.affyexpr;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class AffyExperiment {
    
    public static void main(String[] args) { 
        File probeFile = new File(args[0]);
        try {
            AffyProbeSet probeSet = new AffyProbeSet(probeFile);
            File exptDir = new File(args[1]);
            AffyExperiment[] array = loadExperiments(probeSet, exptDir);

			for(int i = 0; i < array.length; i++) { 
				System.out.println("(" + array[i].size() + ") " + array[i].getBaseFile().getName());
			}

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public static AffyExperiment[] loadExperiments(AffyProbeSet probeSet, File dir) throws IOException { 
        File[] files = dir.listFiles();
        AffyExperiment[] array = new AffyExperiment[files.length];
        for(int i = 0; i < files.length; i++) { 
            array[i] = new AffyExperiment(probeSet, files[i]);
        }
        return array;
    }
    
    private AffyProbeSet probes;
    private Map<AffyProbe,AffyMeasurement> measures;
	private File baseFile;

    public AffyExperiment(AffyProbeSet probeSet, File f) throws IOException {
		baseFile = f;
        probes = probeSet;
        measures = new HashMap<AffyProbe,AffyMeasurement>();
        parseFile(f);
    }
    
    public void addFileToExperiment(File f) throws IOException { 
        parseFile(f);
    }
    
    private void parseFile(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = null;
        
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0 && !line.startsWith("Probe")) { 
                AffyMeasurement m = new AffyMeasurement(probes, line);
                if(m.getProbe() != null) { 
                    if(measures.containsKey(m.getProbe())) { 
                        throw new IllegalArgumentException(line); 
                    }
                    measures.put(m.getProbe(), m);
                }
            }
        }
        br.close();        
    }
   
	public File getBaseFile() { return baseFile; }
    public int size() { return measures.size(); }
    public Collection<AffyProbe> getMeasuredProbes() { return measures.keySet(); }
    public boolean containsProbe(AffyProbe probe) { return measures.containsKey(probe); }
    public AffyMeasurement getMeasurement(AffyProbe ap) { return measures.get(ap); }
    
    public Vector<Double> lookupGeneSymbolValues(String v) { 
        return collectValues(probes.lookupGeneSymbol(v)); 
    }
    public Vector<String> lookupGeneSymbolCalls(String v) { 
        return collectCalls(probes.lookupGeneSymbol(v)); 
    }
    
    public Vector<Double> lookupGeneNameValues(String v) { 
        return collectValues(probes.lookupGeneName(v)); 
    }
    public Vector<String> lookupGeneNameCalls(String v) { 
        return collectCalls(probes.lookupGeneName(v)); 
    }
    
    public Vector<Double> lookupUnigeneValues(String v) { 
        return collectValues(probes.lookupUnigene(v)); 
    }
    public Vector<String> lookupUnigeneCalls(String v) { 
        return collectCalls(probes.lookupUnigene(v)); 
    }
    
    public Vector<Double> lookupLocusIDValues(String v) { 
        return collectValues(probes.lookupLocusID(v)); 
    }
    public Vector<String> lookupLocusIDCalls(String v) { 
        return collectCalls(probes.lookupLocusID(v)); 
    }
    
    private Vector<Double> collectValues(Collection<AffyProbe> pc) { 
        Vector<Double> values = new Vector<Double>();
		if(pc != null) { 
			for(AffyProbe ap : pc) { 
				if(measures.containsKey(ap)) {
					values.add(measures.get(ap).getValue());
				}
			}
		}
        return values;
    }

    private Vector<String> collectCalls(Collection<AffyProbe> pc) { 
        Vector<String> values = new Vector<String>();
		if(pc != null) { 
			for(AffyProbe ap : pc) { 
				if(measures.containsKey(ap)) {
					values.add(measures.get(ap).getCall());
				}
			}
		}
        return values;
    }
}
