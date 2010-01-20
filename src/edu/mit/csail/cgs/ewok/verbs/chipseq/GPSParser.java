package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

/* 
 * Class to parse out put from GPS Chip-Seq binding peak finder
 */
public class GPSParser {

	/**
	 * Parses data in the GPS output format
	 * @param filename name of the file containing the data
	 * @return a List of hit objects
	 */
	public static List<GPSPeak> parseGPSOutput(String filename, Genome g) {
		ArrayList<GPSPeak> results = new ArrayList<GPSPeak>();

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
				GPSPeak hit = GPSParser.parseLine(g, line, 0);
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
	 * 
	 * old format
	 * Position	Rank_Sum	Strength	EM_Posi	Shape	Shape_Z	Shape_Param	ShapeAsymmetry	IpStrength	CtrlStrength	Q_value_log10	MixingProb	NearestGene	Distance	
	 * 8:20401711	47581	200.0	 8:20401711	1.30	1.6223	7.454790	0.33			200.0		4.8				47.53			0.9745		Hoxa1		2147483647	
	   New format
	   Position		IpStrength	Shape	CtrlStrength	Q_value_log10	P_value_log10	UnaryEvent	NearestGene	Distance	Alpha
	   18:75725340	230.5		-0.13	0.4				62.51			67.17			1			NONE		2147483647	9.0	

	 * @param gpsLine a line of text representing a hit
	 * @return a hit object containing the data from the specified line
	 */
	private static GPSPeak parseLine(Genome g, String gpsLine, int lineNumber) {
		GPSPeak peak;
		String[] t = gpsLine.split("\t");
		if (t.length == 10) {
			try { 
				Region r = Region.fromString(g, t[0]);
//				Region em_pos = Region.fromString(g, t[3]);
//				GPSPeak(Genome g, String chr, int pos, int EM_pos, double strength, 
//						double controlStrength, double qvalue, double shape, double shapeZ)
				peak = new GPSPeak(g, r.getChrom(), r.getStart(), 
						Double.parseDouble(t[1]), Double.parseDouble(t[3]), Double.parseDouble(t[4]), 
						Double.parseDouble(t[5]), Double.parseDouble(t[2]), Integer.parseInt(t[6]), t[7], Integer.parseInt(t[8]));
			}
			catch (Exception ex) {
				//logger.error("Parse error on line " + lineNumber + ".", ex);
				return null;
			}
		}
		else {
			//logger.error("Line " + lineNumber + " has " + t.length + " tokens.");
			System.err.println("Line " + lineNumber + " has " + t.length + " tokens.");
			return null;
		}
		return peak;
	}
}
