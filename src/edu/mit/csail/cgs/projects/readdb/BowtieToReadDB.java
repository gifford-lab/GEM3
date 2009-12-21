package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;

/**
 * Reads Bowtie output on stdin.
 * Produces a file on stdout in the format expected by ImportHits.
 * The weight for a hit is 1/(# of hits for that read)
 * 
 * Options:	--nosuboptimal (flag to only take the hits with the minimum number of mismatches)
 * 			--uniquehits (flag to only print 1:1 read to hit mappings)
 * 
 * nosuboptimal is applied before uniquehits
 */

public class BowtieToReadDB {

	private static int readLength=-1;
	
    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("s","nosuboptimal",false,"do not include hits whose score is not equal to the best score for the read");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	boolean uniqueOnly = cl.hasOption("uniquehits");
    	boolean filterSubOpt = cl.hasOption("nosuboptimal");
    	
        ArrayList<String[]> lines = new ArrayList<String[]>();

        String line;
        String lastRead = "";        
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\t");
            if (!pieces[0].equals(lastRead)) {
                printLines(lines, uniqueOnly, filterSubOpt);
                lines.clear();
            }
            lines.add(pieces);
            lastRead = pieces[0];
        }
        printLines(lines, uniqueOnly, filterSubOpt);
    }
    
    /* ImportHits (where the printed lines are eventually used) takes
     * the average between the start and stop positions in the
     * coordinate.  I think it is most correct to have this value
     * reflect the 5' coordinate for each hit. This is why I do the
     * calculation gymnastics for the "coord" String below (even
     * though this will slow imports down a bit) and also why the
     * start and stop positions are the same.
     */ 

    public static void printLines(ArrayList<String[]> lines, boolean uniqueOnly, boolean filterSubOpt) {
        double numHits=0;
        ArrayList<String[]> linesToPrint;
    	//First count the number of (valid) hits
        if(filterSubOpt && lines.size()>1){
    		linesToPrint = new ArrayList<String[]>();
    		int minMis=Integer.MAX_VALUE;
    		//Find minimum number of mismatchs
    		int [] mismatches = new int [lines.size()];
    		int i=0;
    		for(String[] pieces : lines){
    			int mis=0;
    			if(pieces.length>7 && pieces[7].length()>1)
    				mis = pieces[7].split(",").length;
    			mismatches[i]=mis;
    			if(mis<minMis)
    				minMis=mis;
    			i++;
    		}
    		//Only add minimum mismatch hits
    		i=0;
    		for(String[] pieces : lines){
    			if(mismatches[i]==minMis)
    				linesToPrint.add(pieces);
    			i++;
    		}
    		numHits = linesToPrint.size();
    	}else{
    		linesToPrint=lines;
    		numHits=lines.size();
    	}
    	
    	//Now add the hits!
    	if(!uniqueOnly || linesToPrint.size()==1){
	    	double weight = 1.0 / numHits;
	        for (String[] pieces : linesToPrint) {
	        	if(readLength==-1)
	        		readLength=pieces[4].length();
	        	
	        	String coord = pieces[1].equals("+") ? pieces[3] : Integer.toString((Integer.valueOf(pieces[3])+readLength-1)); 
	            System.out.println(String.format("%s\t%s\t%s\t%f",
	                                             pieces[2],
	                                             coord,
	                                             pieces[1],
	                                             weight));
	        }
    	}
    }
}
