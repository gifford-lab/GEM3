package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPSFastaWriter{  
	// --species "Homo sapiens;hg19" --window 200 --top -1 --no_cache  --root C:\Data\workspace\gse [--expts expt_done.txt] --write_read_coverage
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
    
     SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	seqgen.useCache(!flags.contains("no_cache"));
	boolean skip_repeat = !flags.contains("allow_repeat");
	boolean write_read_coverage = flags.contains("write_read_coverage");
	
	List<String> names = new ArrayList<String>();
	File dir = new File(Args.parseString(args, "root", null));
    if (!dir.exists()){
      System.err.println("Please provide root of GEM/GPS analysis folders, '--root root_path' ");
      System.exit(0);
    }	
    
    String expts_file = Args.parseString(args, "expts", null);
    if (expts_file!=null){		// read expt names from a file
    	File listFile = new File(dir, expts_file);
	    if (!listFile.exists()){
	    	System.err.println("Expt list file not exist: "+listFile.getAbsolutePath());
	    	System.exit(0);
	    }
		BufferedReader bin = null;
		try{
	        bin = new BufferedReader(new InputStreamReader(new FileInputStream(listFile)));
	        String line;
	        while((line = bin.readLine()) != null) 
	            names.add(line.trim());
		}
		catch (IOException e){
			System.err.println("Error in reading expt list file, "+listFile.getAbsolutePath());
			e.printStackTrace(System.err);
		}
    }
    else{
		File[] children = dir.listFiles();
		for (int i=0;i<children.length;i++){
			File child = children[i];
			if (child.isDirectory())
				names.add(child.getName());
		}
    }
	
    // load GPS results
	for (String exptName : names){
		System.out.println("Writing sequences and related info from events for "+exptName);
	    File gpsFile = new File(new File(exptName), exptName+"_1_GPS_significant.txt");
	    if (!gpsFile.exists()){
	    	System.err.println("GPS file not exist: "+gpsFile.getAbsolutePath());
	        continue;
	    }
	    List<GPSPeak> gpsPeaks =null;
	    try{
	    	gpsPeaks= GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
	    }
	    catch(IOException e){
	    	System.err.println("Error reading/parsing GPS file: "+gpsFile.getAbsolutePath());
	        continue;
	    }
	    
	    // use exptName.report.txt file to get the data strings to read db
	    DeepSeqExpt chipSeq = null;
	    if (write_read_coverage){
	    	ArrayList<String> readDbStrings = new ArrayList<String>();
		    File reportFile = new File(new File(exptName), exptName+".report.txt");
		    if (!reportFile.exists()){
		    	System.err.println("Report file not exist: "+reportFile.getAbsolutePath());
		        continue;
		    }
			BufferedReader bin = null;
			try{
		        bin = new BufferedReader(new InputStreamReader(new FileInputStream(reportFile)));
		        String line;
		        while((line = bin.readLine()) != null) { 
		            String[] f=line.split("--rdbexptE1 \"");
		            if (f.length>1){
		            	for (int i=1;i<f.length;i++){
		            		String field = f[i];
		            		int end = field.indexOf("\"");
		            		readDbStrings.add(field.substring(0,end));
		            	}
		            	break;
		            }
		        }
			}
			catch (IOException e){
				System.err.println("Error in reading report file, "+reportFile.getAbsolutePath());
				e.printStackTrace(System.err);
			}
			
		    // get the read data
		    ArrayList<ChipSeqLocator> locators = new ArrayList<ChipSeqLocator>();
	        for (String dbStr: readDbStrings) {
	            String[] pieces = dbStr.split(";");
	            if (pieces.length == 2) {
	            	locators.add(new ChipSeqLocator(pieces[0], pieces[1]));
	            } else if (pieces.length == 3) {
	            	locators.add(new ChipSeqLocator(pieces[0], pieces[1], pieces[2]));
	            } else {
	                throw new RuntimeException("Couldn't parse a ChipSeqLocator from " + dbStr);
	            }
	        }
	        if (locators.isEmpty())
	        	continue;
	        chipSeq = new DeepSeqExpt(genome, locators, "readdb", -1);
	    }        
	    
	    if (top==-1)				// if top = -1, use all peaks
	    	top = gpsPeaks.size();
	    else
	    	top = Math.min(top, gpsPeaks.size());
	    int count=1;
	    StringBuilder sb = new StringBuilder();
	    StringBuilder hms_summit_sb = new StringBuilder();
	    StringBuilder hms_readcoverage_sb = new StringBuilder();
	    StringBuilder chipmunk_sb = new StringBuilder();
	    eachPeak: for (GPSPeak p: gpsPeaks){
	    	if (count>top)
	    		break;
	    	int start = p.getLocation()-window/2;
	    	if (start<0)
	    		start=0;
	    	int end = start+window-1;
	    	if (end>=genome.getChromLength(p.getChrom()))
	    		continue;
	    	Region r = new Region(genome, p.getChrom(), start, end);
	    	String seq = seqgen.execute(r);
	    	if (skip_repeat){
				for (char c:seq.toCharArray())
					if (Character.isLowerCase(c) || c=='N')
						continue eachPeak;
	    	}
	    	//	passed the repeat check, output
	    	sb.append(">seq_").append(count).append(" ").append(exptName).append(" ").append(r.toString()).append("\n");
	    	sb.append(seq).append("\n");
	    	hms_summit_sb.append(window/2).append("\n");
	    	if (write_read_coverage){
	    		// HMS coverage is based on read itself only
	    		double[] coverage = new double[r.getWidth()]; 
	    		int origin = r.getStart();
	    		List<ReadHit> hits = chipSeq.loadHits(r);
	    		for (ReadHit h : hits){
	    			for (int i=h.getStart();i<=h.getEnd();i++){
	    				int idx = i-origin;
	    				if (idx>=0 && idx<r.getWidth())
	    					coverage[idx]++;
	    			}
	    		}
	    		hms_readcoverage_sb.append(">seq_").append(count).append(" ").append(exptName).append(" ").append(r.toString()).append("\n");
	    		for (int i=0;i<coverage.length;i++){
	    			hms_readcoverage_sb.append(String.format("(%d) %.2f ", i, coverage[i]));
	    		}
	    		hms_readcoverage_sb.append("\n");
	    		
	    		// ChIPMunk coverage is based on extended reads (here extend to 200bp)
	    		coverage = new double[r.getWidth()];
	    		int readExtendedLength = 200;
	    		for (ReadHit h : hits){
	    			if (h.getStrand()=='+'){
		    			for (int i=h.getStart();i<h.getStart()+readExtendedLength;i++){
		    				int idx = i-origin;
		    				if (idx>=0 && idx<r.getWidth())
		    					coverage[idx]++;
		    			}
	    			}
	    			else{
		    			for (int i=h.getEnd()-readExtendedLength+1;i<=h.getEnd();i++){
		    				int idx = i-origin;
		    				if (idx>=0 && idx<r.getWidth())
		    					coverage[idx]++;
		    			}
	    			}
	    		}
	    		chipmunk_sb.append("> ");
	    		for (int i=0;i<coverage.length;i++){
	    			chipmunk_sb.append(String.format("%.2f ", coverage[i]));
	    		}
	    		chipmunk_sb.append("\n").append(seq).append("\n");
	    	}
	    	count++;
	    }
	    if (write_read_coverage){
		    CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.summit.txt", hms_summit_sb.toString());
		    CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.fasta", sb.toString());
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.basecover.txt", hms_readcoverage_sb.toString());
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_ChIPMunk.peak.txt", chipmunk_sb.toString());
	    }
	    else
	    	CommonUtils.writeFile(exptName+"_"+window+"bp.fasta", sb.toString());
	    System.out.println(exptName+" is processed, "+(count-1)+" sequences has been written.");
	  }
  }
}
