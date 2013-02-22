package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

/* 
 * Class to parse out put from GPS Chip-Seq binding peak finder
 */
public class GPSParser {
	
	public static void main(String[] args) throws IOException {
		Genome genome=null;
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
	      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
	    String filename = Args.parseString(args,"GPS",null); 
	    List<GPSPeak> peaks = parseGPSOutput(filename, genome);
	    
	    String outFormat = Args.parseString(args,"of",null); 
	    if (outFormat==null){
		    System.out.println(peaks.size());
		    System.out.println(GPSPeak.toGPS_Header());
		    System.out.println(peaks.get(0).toGPS());
		    System.out.println(GPSPeak.toGPS_short_Header());
		    System.out.println(peaks.get(0).toGPS_short());
	    }
	    if (outFormat.equals("narrow_peak")){
	    	int count=0;
	    	for (int i=0;i<peaks.size();i++){
	    		GPSPeak p = peaks.get(i);
	    		count++;
	    		if (p.getQvalue()<999)
	    			break;
	    	}
//	    	for (int i=0;i<100;i++){
	    	for (int i=0;i<peaks.size();i++){
	    		GPSPeak p = peaks.get(i);
	    		double score = p.getQvalue()<999?p.getQvalue():(1000+count-i);
	    		System.out.println(p.toNarrowPeak(score));
	    	}
	    }
	}
	/**
	 * Parses data in the GPS output format
	 * @param filename name of the file containing the data
	 * @return a List of hit objects
	 */
	public static List<GPSPeak> parseGPSOutput(String filename, Genome g) throws IOException {
        return parseGPSOutput(new FileInputStream(filename), g);
    }


	public static List<GPSPeak> parseGPSOutput(InputStream stream, Genome g) throws IOException {
		ArrayList<GPSPeak> results = new ArrayList<GPSPeak>();

		BufferedReader bin = null;
		int count = 0;
        bin = new BufferedReader(new InputStreamReader(stream));
		
        String line;
        while((line = bin.readLine()) != null) { 
            line = line.trim();
            String[] f=line.split("\t");
            if (f[0].equals("Position")||f[0].equals("chr")){
                continue;
            }
            try {
                GPSPeak hit = GPSParser.parseLine(g, line, 0);
                if (hit!=null) {
                    results.add(hit);
                }
            } catch (RuntimeException e) {
                System.err.println("Parsing line " + line);
                System.err.println(e.toString());
                throw(e);
            }				
            count++;
        }			
        if (bin != null) {
            bin.close();
        }
		return results;
	}

	/**
	 * Parse a single line of text into a GPSPeak object
	 * New format
	   Position	      IP	Control	   Fold	Q_-lg10	P_-lg10	IPvsEMP	IPvsCTR	Joint	NearestGene	Distance	Alpha	EM_Position
	   4:117003586	  496.8	    7.4	   67.1	129.052	135.017	 -0.772	  4.490	1	NONE	2147483647	    7.8	4:117003586
	 * @param gpsLine a line of text representing a hit
	 * @return a hit object containing the data from the specified line
	 */
	public static GPSPeak parseLine(Genome g, String gpsLine, int lineNumber) {
		GPSPeak peak;
		String[] t = gpsLine.split("\t");
		if (t.length == 14 || t.length == 15 ) {		// with kmer info
	    	// GPS output format 2011-07-25	
	    	// Position	     IP	Control	   Fold	Expectd	Q_-lg10	P_-lg10	P_poiss	IPvsEMP	IPvsCTR	Kmer	Count	Strength	BoundSequence	EnrichedHGP
	    	Region r = Region.fromString(g, t[0]);
            peak = new GPSPeak(g, r.getChrom(), r.getStart(), 
                    Double.parseDouble(t[1]), Double.parseDouble(t[2]), Double.parseDouble(t[4]), Double.parseDouble(t[5]), 
                    Math.pow(10,-1*Double.parseDouble(t[6])), Double.parseDouble(t[6]), Double.parseDouble(t[7]), Double.parseDouble(t[8]), Double.parseDouble(t[9]),
                    t[10], (int)Double.parseDouble(t[11]), t[12].charAt(0), t[13]);
	    } else if (t.length == 13 ) {		// with kmer info
	    	// GPS output format 2012-03	
	    	// Position	     IP	Control	   Fold	Expectd	Q_-lg10	P_-lg10	P_poiss	IPvsEMP	IPvsCTR	Kmer	Count	Strength	BoundSequence
	    	Region r = Region.fromString(g, t[0]);
            peak = new GPSPeak(g, r.getChrom(), r.getStart(), 
                    Double.parseDouble(t[1]), Double.parseDouble(t[2]), Double.parseDouble(t[4]), Double.parseDouble(t[5]), 
                    Math.pow(10,-1*Double.parseDouble(t[6])), Double.parseDouble(t[6]), Double.parseDouble(t[7]), Double.parseDouble(t[8]), Double.parseDouble(t[9]),
                    t[10], (int)Double.parseDouble(t[11]), t[12].charAt(0), "");
	    } else {
            throw new RuntimeException("Invalid number of fields (" + t.length + ") on line " + lineNumber + ": " + gpsLine);
        }
        return peak;
    }
}
