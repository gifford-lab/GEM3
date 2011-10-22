package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPSFastaWriter{  
	
  public static void main(String[] args){
    ArgParser ap = new ArgParser(args);
    Set<String> flags = Args.parseFlags(args);
    Genome genome=null;
    
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
      }  
    } catch (NotFoundException e) {
      e.printStackTrace();
    }
    
    int window = Args.parseInteger(args, "window", 100);
    int top = Args.parseInteger(args, "top", 500);
    
    // load GPS results
    String GPSfileName = Args.parseString(args, "GPS", null);
    if (GPSfileName==null){
      System.err.println("Please provide GPS file, '--GPS file' ");
      System.exit(0);
    }
    File gpsFile = new File(GPSfileName);
    if (!gpsFile.exists()){
        System.err.println("Can not find GPS file : "+GPSfileName);
        System.exit(0);
    }

    SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	seqgen.useCache(!flags.contains("no_cache"));

	GPSfileName = gpsFile.getName();
    String exptName = GPSfileName.substring(0,GPSfileName.indexOf("_1_GPS_significant.txt"));
    List<GPSPeak> gpsPeaks =null;
    try{
    	gpsPeaks= GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
    }
    catch(IOException e){
    	System.err.println("Error reading/parsing GPS file: "+GPSfileName);
        System.exit(0);
    }
    top = Math.min(top, gpsPeaks.size());
    int count=0;
    StringBuilder sb = new StringBuilder();
    for (GPSPeak p: gpsPeaks){
    	if (count>top)
    		break;
    	int start = p.getLocation()-window/2;
    	if (start<0)
    		start=0;
    	int end = start+window-1;
    	if (end>=genome.getChromLength(p.getChrom()))
    		continue;
    	Region r = new Region(genome, p.getChrom(), start, end);
    	sb.append(">").append(exptName).append("\t#").append(count).append("\t").append(r.toString()).append("\n");
    	sb.append(seqgen.execute(r)).append("\n");
    	count++;
    }
    CommonUtils.writeFile(exptName+"_"+window+"bp.fasta", sb.toString());
  }
  
}
