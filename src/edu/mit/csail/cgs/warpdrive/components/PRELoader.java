/*
 * Created on Oct 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive.components;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PRELoader {
    
    private int width;
    private Vector<MotifScanner> scanners;
    private Map<String,Double> singleWeights;
    private Map<String,Map<String,Double>> pairWeights;
    
    public PRELoader(int w, File f) throws IOException { 
        width = w;
        singleWeights = new HashMap<String,Double>();
        pairWeights = new HashMap<String,Map<String,Double>>();
        scanners = new Vector<MotifScanner>();

        BufferedReader br = new BufferedReader(new FileReader(f));
        
        boolean readingMotifs = true;
        Pattern mp = Pattern.compile("motif\\s+([^\\s]+)\\s+([^\\s]+)\\s*");
        Pattern swp = Pattern.compile("weight\\s+([^\\s:]+)\\s+([^\\s]+)\\s*");
        Pattern pwp = Pattern.compile("weight\\s+([^\\s]+):([^\\s]+)\\s+([^\\s]+)\\s*");
            
        String line;
        while((line = br.readLine()) != null) { 
            if((line = line.trim()).length() > 0) { 
                if(readingMotifs) { 
                    if(line.startsWith("---")) {
                        readingMotifs = false;
                        
                        for(MotifScanner ms : scanners) { 
                            singleWeights.put(ms.getName(), 0.0);
                            pairWeights.put(ms.getName(), new HashMap<String,Double>());
                            for(MotifScanner ms2 : scanners) { 
                                pairWeights.get(ms.getName()).put(ms2.getName(), 0.0);
                            }
                        }
                        
                    } else { 
                        Matcher m = mp.matcher(line);
                        if(m.matches()) { 
                            String name = m.group(1);
                            String regex = m.group(2);
                            RegexMotifScanner scanner = new RegexMotifScanner(name, regex);
                            scanners.add(scanner);
                        }
                    }
                } else { 
                    Matcher m = swp.matcher(line);
                    if(m.matches()) { 
                        String ms = m.group(1);
                        double weight = Double.parseDouble(m.group(2));
                        singleWeights.put(ms, weight);
                        System.out.println(String.format("%s --> %f", ms, weight));
                    } else { 
                        m = pwp.matcher(line);
                        if(m.matches()) { 
                            String ms1 = m.group(1);
                            String ms2 = m.group(2);
                            double weight = Double.parseDouble(m.group(3));
                            pairWeights.get(ms1).put(ms2,weight);
                            pairWeights.get(ms2).put(ms1,weight);
                            System.out.println(String.format("%s,%s --> %f", ms1, ms2, weight));
                        }
                    }
                }
            } 
        }
        br.close();
    }
    
    public PREScanner getScanner() { 
        return new PREScanner(width, scanners, singleWeights, pairWeights);
    }
}
