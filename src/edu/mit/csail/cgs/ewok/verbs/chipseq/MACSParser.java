package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.species.Genome;


public class MACSParser {

	/**
	 * Parses data in the MACS output, e.g.
	 * # d = 124
	 * chr	start	end	length	summit	tags	-10*log10(pvalue)	fold_enrichment	FDR(%)
	 * chr1	4481542	4484063	2522	276	272	2629.53	47.51	0.00

	 * @param filename name of the file containing the data
	 * @return a List of hit objects
	 */
	public static List<MACSPeakRegion> parseMACSOutput(String filename, Genome g) {
		ArrayList<MACSPeakRegion> results = new ArrayList<MACSPeakRegion>();

		FileReader in = null;
		BufferedReader bin = null;
		try {
			in = new FileReader(filename);
			bin = new BufferedReader(in);

			
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
	            String[] f=line.split("\t");
	            if (line.charAt(0)=='#'||f[0].equals("chr")){
	            	continue;
	            }
				MACSPeakRegion hit = MACSParser.parseLine(g, line, 0);
				if (hit!=null)
					results.add(hit);
			}
		}
		catch(IOException ioex) {
			//logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				//logger.error("Error closing buffered reader", ioex2);
			}			
		}
		return results;
	}



	/**
	 * Parse a single line of text into a hit object
	 * @param macsLine a line of text representing a hit
	 * @return a hit object containing the data from the specified line
	 */
	private static MACSPeakRegion parseLine(Genome g, String macsLine, int lineNumber) {
		MACSPeakRegion macs;
		String[] t = macsLine.split("\t");
		if (t.length == 9) {
			try { 
				macs = new MACSPeakRegion(g, t[0], 
						Integer.parseInt(t[1]),
						Integer.parseInt(t[2]),
						Integer.parseInt(t[4]),
						Integer.parseInt(t[5]),
						Double.parseDouble(t[6]),
						Double.parseDouble(t[7]),
						Double.parseDouble(t[8]));
			}
			catch (Exception ex) {
				//logger.error("Parse error on line " + lineNumber + ".", ex);
				return null;
			}
		}
		else {
			//logger.error("Line " + lineNumber + " has " + tokens.length + " tokens.");
			return null;
		}
		return macs;
	}
}
