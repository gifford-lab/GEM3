package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class GPSFastaWriter{  
	// --species "Homo sapiens;hg19" --expts expt_done.txt [--root C:\Data\workspace\gse] --window 200 --top -1 --no_cache  
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
    
    boolean wantMEME = flags.contains("meme");
    boolean wantGEM = flags.contains("gem");
    boolean wantHMS = flags.contains("hms");
    boolean wantChIPMunk = flags.contains("chipmunk");
    
    int window = Args.parseInteger(args, "window", 100);
    int top = Args.parseInteger(args, "top", 500);
    int k_neg_dist = Args.parseInteger(args, "k_neg_dist", 300);
    
    SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	seqgen.useCache(!flags.contains("no_cache"));
	seqgen.useLocalFiles(!flags.contains("use_db_genome"));
	
	boolean skip_repeat = flags.contains("no_repeat");
	
	List<String> names = new ArrayList<String>();
	File root = new File(Args.parseString(args, "root", null));
    if (!root.exists()){
      System.err.println("Please provide root of GEM/GPS analysis folders, '--root root_path' ");
      System.exit(0);
    }	
    
    String expts_file = Args.parseString(args, "expts", null);
    if (expts_file!=null){		// read expt names from a file
    	File listFile = new File(expts_file);
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
    else{		// if no expt name file, try all sub-folders
		File[] children = root.listFiles();
		for (int i=0;i<children.length;i++){
			File child = children[i];
			if (child.isDirectory())
				names.add(child.getName());
		}
    }
	
    // load GPS results
	for (String exptName : names){
		System.out.print("Processing sequences and related info for "+exptName+", ");
		File folder = new File(new File(root, exptName), exptName+"_outputs");
	    File gpsFile = new File(folder, exptName+"_1_GEM_events.txt");
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
	    if (wantHMS || wantChIPMunk){
	    	ArrayList<String> readDbStrings = new ArrayList<String>();
		    File reportFile = new File(new File(root, exptName), exptName+".report.txt");
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
	    } // prepare read DB access        
	    
	    int seqNum = 0;
	    if (top==-1)				// if top = -1, use all peaks
	    	seqNum = gpsPeaks.size();
	    else
	    	seqNum = Math.min(top, gpsPeaks.size());
	    int count=1;
	    StringBuilder sb = new StringBuilder();
	    StringBuilder gem_sb = new StringBuilder();
	    StringBuilder gem_neg_sb = new StringBuilder();
	    StringBuilder hms_summit_sb = new StringBuilder();
	    StringBuilder hms_readcoverage_sb = new StringBuilder();
	    StringBuilder chipmunk_sb = new StringBuilder();
	    ArrayList<Region> negativeRegions = new ArrayList<Region>();

	    ArrayList<Region> expandedRegions = new ArrayList<Region>(); // regions expanded from positive positions
	    eachPeak: for (GPSPeak p: gpsPeaks){
	    	if (count>seqNum)
	    		break;
	    	int start = p.getLocation()-window/2;
	    	if (start<0)
	    		start=0;
	    	int end = start+window-1;
	    	if (end>=genome.getChromLength(p.getChrom()))
	    		continue;
	    	Region r = new Region(genome, p.getChrom(), start, end);
	    	if (wantGEM){
	    		r = p.expand(window/2);
	    		if (r.getWidth()!=2*(window/2)+1)		// if at the end of chromosome, skip
	    			continue;
	    	}
	    	String seq = seqgen.execute(r);
	    	// negative region for GEM
			// In proximal regions, but excluding binding regions
			String chr = p.getChrom();
			int chrLength = genome.getChromLength(chr)-1;
			start = p.getLocation()+ k_neg_dist;
			if ( start+window>=chrLength)
				continue;
			Region r_neg = new Region(genome, chr, start, start+window);
			negativeRegions.add(r_neg);

			for (char c:seq.toCharArray())
				if ((Character.isLowerCase(c) && skip_repeat) || c=='N')
					continue eachPeak;
			
	    	//	passed the repeat check, output
	    	sb.append(">seq_").append(count).append(" ").append(exptName).append(" ").append(r.toString()).append("\n");
	    	sb.append(seq).append("\n");
	    	
	    	gem_sb.append(String.format(">seq_%d %.1f %s %s\n", count, p.getStrength(), exptName, r.toString()));
	    	gem_sb.append(seq).append("\n");

	    	hms_summit_sb.append(window/2).append("\n");
	    	
	    	List<ReadHit> hits = null;
	    	if (wantHMS){
	    		// HMS coverage is based on read itself only
	    		double[] coverage = new double[r.getWidth()]; 
	    		int origin = r.getStart();
	    		hits = chipSeq.loadHits(r);
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
	    		
	    	}
	    	if (wantChIPMunk){
	    		// ChIPMunk coverage is based on extended reads (here extend to 200bp)
	    		double[] coverage = new double[r.getWidth()];
	    		int origin = r.getStart();
	    		int readExtendedLength = 200;
	    		if (hits==null)
	    			hits = chipSeq.loadHits(r);
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
	    
	    // GEM negative regions, excluding positive regions

		negativeRegions = Region.filterOverlapRegions(negativeRegions, expandedRegions);
		int neg_count=0;
		for (Region r: negativeRegions){			
	    	String neg_seq = seqgen.execute(r);
	    	gem_neg_sb.append(String.format(">neg_seq_%d %s %s\n", neg_count, exptName, r.toString()));
	    	gem_neg_sb.append(neg_seq).append("\n");
	    	neg_count++;
		}
		
	    if (wantHMS || wantChIPMunk){
			chipSeq.closeLoaders();
			chipSeq=null;
			System.gc();
		    CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.summit.txt", hms_summit_sb.toString());
		    CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.fasta", sb.toString());
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_HMS.basecover.txt", hms_readcoverage_sb.toString());
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_ChIPMunk.peak.txt", chipmunk_sb.toString());
	    }
	    if (wantMEME)
	    	CommonUtils.writeFile(exptName+"_"+window+"bp.fasta", sb.toString());
	    if (wantGEM){
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_GEM.fasta", gem_sb.toString());
	    	CommonUtils.writeFile(exptName+"_"+window+"bp_GEM_neg.fasta", gem_neg_sb.toString());
	    }
	    
	    System.out.println((count-1)+" sequences.");
	  }
  }
}
