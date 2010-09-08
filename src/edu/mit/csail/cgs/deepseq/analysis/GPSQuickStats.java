package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPSQuickStats{
  
  private Genome genome;
  private Organism org;
  private String[] args;
  private WeightMatrix motif = null;
  private double motifThreshold;
  private int motif_window;
  private String GPSPrefix;
  
  public static void main(String[] args) throws IOException {
    
	GPSQuickStats analysis = new GPSQuickStats(args);
	analysis.execute();
  }
  
  public GPSQuickStats(String[] args) {
    this.args = args;
    ArgParser ap = new ArgParser(args);
    
    try {
      Pair<Organism, Genome> pair = Args.parseGenome(args);
      if(pair==null){
        //Make fake genome... chr lengths provided???
        if(ap.hasKey("geninfo")){
          genome = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
            }else{
              System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");;System.exit(1);
            }
      }else{
        genome = pair.cdr();
        org = pair.car();
      }  
    } catch (NotFoundException e) {
      e.printStackTrace();
    }
	// load motif
    String motifString = Args.parseString(args, "motif", null);
    if (motifString!=null){
		Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
		motif = wm.car();
		motifThreshold = wm.cdr();
    }
    motif_window = Args.parseInteger(args, "windowSize", 100);
  }
  
  private void execute() throws IOException {
	GPSPrefix = Args.parseString(args, "prefix", null);
	if (GPSPrefix==null){
	  System.err.println("GPS prefix not found!");
	  System.exit(0);
	}
	
	StringBuilder sb= new StringBuilder();
	sb.append(runOutputAnalysis(GPSPrefix+"_GPS_significant.txt"))
	.append(runOutputAnalysis(GPSPrefix+"_GPS_insignificant.txt"))
	.append(runOutputAnalysis(GPSPrefix+"_GPS_filtered.txt"));
	CommonUtils.writeFile("Analysis_"+GPSPrefix+"_Stats.txt", sb.toString());
  }
  
  private String runOutputAnalysis(String GPSfileName) throws IOException{  
	long tic = System.currentTimeMillis();
    File gpsFile = new File(GPSfileName);
    List<GPSPeak> gpsPeaks = null;
    try{
    	gpsPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
    }
    catch (IOException e){
    	return "";
    }
//	Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
//	    public int compare(GPSPeak o1, GPSPeak o2) {
//	        return o1.compareByPValue(o2);
//	    }
//	});
	
	GPSOutputAnalysis gpsAnalysis = new GPSOutputAnalysis(genome, motif, motifThreshold, gpsPeaks, 
			"Analysis_"+GPSfileName.substring(0, GPSfileName.length()-4), motif_window);
	String stats = gpsAnalysis.printMotifHitList();
	String msg = "\n"+CommonUtils.timeElapsed(tic)+ "\n";
	System.out.print(msg);
	return stats+msg;
  }
  
  
  
}
 