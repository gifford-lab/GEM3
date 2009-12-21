package edu.mit.csail.cgs.utils.parsing.affyexpr;

import affymetrix.fusion.cdf.*;
import affymetrix.fusion.cel.*;
import java.io.*;
import java.util.*;

public class AffyFusionCelReader {
	
	public static void main(String[] args) { 
		String cdfFile = args[0];
        String celFile = args[1];
		AffyFusionCelReader reader = new AffyFusionCelReader(cdfFile, celFile);
        
        if(reader.read()) { 
            Map<String,Double> probeValues = reader.translate();
			System.out.println("# Values: " + probeValues.size());
        } else { 
            System.err.println("Couldn't read CDF or CEL file.");
        }
	}
	
	private FusionCELData cel;
	private FusionCDFData cdf;

	public AffyFusionCelReader(String cdfFilename, String celFilename) {
		cel= new FusionCELData();
		cel.setFileName(celFilename);
		cdf = new FusionCDFData();
		cdf.setFileName(cdfFilename);
        
        if(!cdf.exists() || !cel.exists()) { throw new IllegalArgumentException("CDF or CEL file doesn't exist."); }
	}
    
    public boolean read() { 
        return cdf.read() && cel.read();
    }
    
    public Map<String,String> getParameters() { 
        Map<String,String> params = new HashMap<String,String>();

        Vector paramVector = cel.getParameters();
        for(Object pv : paramVector) { 
            
        }
        
        return params;
    }
    
    public Map<String,Double> translate() {
        Map<String,Double> probeMap = new HashMap<String,Double>();
        
        FusionCDFProbeSetInformation probeSetInfo = new FusionCDFProbeSetInformation();
        FusionCDFProbeGroupInformation groupInfo = new FusionCDFProbeGroupInformation();
        FusionCDFProbeInformation probeInfo = new FusionCDFProbeInformation();

        int numProbeSets = cdf.getHeader().getNumProbeSets();
        for(int probeSet = 0; probeSet < numProbeSets; probeSet++) { 
            
            cdf.getProbeSetInformation(probeSet, probeSetInfo);
            
            int numGroups = probeSetInfo.getNumGroups();
            String probeSetName = cdf.getProbeSetName(probeSet);
            //System.out.println(probeSetName);
            
            double sum = 0.0;
            int n = 0;
            
            for(int group = 0; group < numGroups; group++) { 
                probeSetInfo.getGroup(group, groupInfo);
                
                int numCells = groupInfo.getNumCells();
                
                for(int cell = 0; cell < numCells; cell++) { 
                    groupInfo.getCell(cell, probeInfo);
                    double int1 = cel.getIntensity(probeInfo.getX(), probeInfo.getY());
                    //System.out.println("\t" + group + "," + cell + " : " + int1);
                    
                    n += 1;
                    sum += int1;
                }
                
            }
            
            if(probeMap.containsKey(probeSetName)) { 
                System.err.println("Duplicate Name: " + probeSetName); 
            } else { 
                probeMap.put(probeSetName, n > 0 ? sum / (double)n : 0.0);
            }
            
            //System.out.println();
        }
        
        return probeMap;
    }
	

}
