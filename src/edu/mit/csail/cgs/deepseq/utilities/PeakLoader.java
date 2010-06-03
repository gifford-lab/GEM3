package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.tools.utils.Args;


/**
 * This class loads the peaks from a StatPeak or a MACS file
 * @author geopapa
 *
 */
public class PeakLoader {
	
	private Genome gen;
	
	/**
	 * Files (one for each different condition-experiment)
	 */
	private File[] files;
	
	/**
	 * Names of the conditions (experiments)
	 */
	private String[] condNames;
	
	/**
	 * map from chromosome names to IDs
	 */
	private Map<String, Integer> chrom2ID=new HashMap<String,Integer>();
	
	/**
	 * map from IDs to chromosome names
	 */
	private Map<Integer,String> id2chrom=new HashMap<Integer,String>();

	/**
	 * map from condition (experiment) names to IDs
	 */
	private Map<String, Integer> condName2ID = new HashMap<String, Integer>();
	
	/**
	 * map from IDs to condition (experiment) names
	 */
	private Map<Integer, String> id2condName = new HashMap<Integer, String>();
	
	/**
	 * format of the files <br>
	 * Supported formats: <tt>StatPeak, MACS</tt> 
	 */
	private String fileFormat;
	
	/**
	 * The peaks corresponding to the enriched regions (as lists)
	 */
	private List<Integer>[][] peakLists;
	
	/**
	 * The peaks corresponding to the enriched regions
	 */
	private int[][][] peaks;
	
	/**
	 * The number of peaks for each condition (experiment)
	 */
	private int[] numPeaks;
	
	
	public PeakLoader(Genome g, String[] filesArg, String fileFormatArg) {
		this(g, filesArg, null, fileFormatArg);
	}
	
	public PeakLoader(Genome g, String[] filesArg, String[] condNamesArg, String fileFormatArg) {
		gen       = g;
		files     = new File[filesArg.length];
		for(int i = 0; i < files.length; i++) { files[i] = new File(filesArg[i]); }
		condNames = condNamesArg;
		fileFormat = fileFormatArg;
		
		initialize();
	}
	
	public PeakLoader(Genome g, List<File> filesArg, String fileFormatArg) {
		this(g, filesArg.toArray(new File[0]), fileFormatArg);
	}

	public PeakLoader(Genome g, List<File> filesArg, String[] condNamesArg, String fileFormatArg) {
		this(g, filesArg.toArray(new File[0]), condNamesArg, fileFormatArg);
	}
	
	public PeakLoader(Genome g, File[] filesArg, String fileFormatArg) {
		this(g, filesArg, null, fileFormatArg);
		
	}
	
	public PeakLoader(Genome g, File[] filesArg, String[] condNamesArg, String fileFormatArg) {
		gen       = g;
		files     = filesArg;
		condNames = condNamesArg;
		fileFormat = fileFormatArg;
				
		initialize();
	}
	
	private void initialize() {
		
		if(condNames == null || condNames.length == 0) {
			condNames = new String[files.length];  for(int i = 0; i < condNames.length; i++) { condNames[i] = files[i].getName(); }
		}
		
		if(condNames.length != files.length)
			throw new IllegalArgumentException("The number of condition (experiment) names" +
					                           "should equal the number of condition (experiment) files.");

		//Initialize the chromosome name lookup tables
		int i=0; 
		for(String c:gen.getChromList()){ chrom2ID.put(c, i); id2chrom.put(i++, c); }
		
		//Initialize the condition (experiment) name lookup tables
		i=0; 
		for(String c:condNames){ condName2ID.put(c, i); id2condName.put(i++, c); }

		peakLists = new ArrayList[condNames.length][gen.getChromList().size()];   for(i = 0; i < peakLists.length; i++) { for(int j = 0; j < peakLists[i].length; j++) { peakLists[i][j] = new ArrayList<Integer>(); } }
		peaks = new int[condNames.length][gen.getChromList().size()][];
		numPeaks = new int[condNames.length];
	}//end of initialize method
	
	
	/**
	 * It loads the peaks of the enriched regions coming either from a StatPeak
	 * or MACS file
	 */
	public void loadPeaks() {
		
		BufferedReader br = null;
		try {
			
			if( fileFormat.equalsIgnoreCase("StatPeak")) {
				int i = 0;
				for(File f:files) {
					
					br = new BufferedReader(new FileReader(f));
					String line;
					while( (line = br.readLine()) != null) {
						if( !line.matches("^[\\w\\d]+:\\d+-\\d+\\s.*") ) { continue; }
						
						String[] words = line.split("\\s");
						String[] peakInfo = words[2].split(":");
						String chromPeak = peakInfo[0]; int posPeak   = Integer.parseInt(peakInfo[1]);
						peakLists[condName2ID.get(condNames[i])][chrom2ID.get(chromPeak)].add(posPeak);
						numPeaks[condName2ID.get(condNames[i])]++;
					}
					i++;
				}//end of for loop
				
			}//end of if condition
			
			else if( fileFormat.equalsIgnoreCase("MACS")) {
				int i = 0;
				for(File f:files) {
					br = new BufferedReader(new FileReader(f));
					String line;
					while( (line = br.readLine()) != null) {
						if( !line.matches("^[Cc][Hh][Rr]\\w*[\\w\\d]+\\s.*") ) { continue; }
						
						String[] words = line.split("\\s");
						
						String chromPeak = words[0]; chromPeak = chromPeak.replaceFirst("chr", "");
						int start  = Integer.parseInt(words[1]); int summit = Integer.parseInt(words[4]);
						int posPeak = start + summit;
						peakLists[condName2ID.get(condNames[i])][chrom2ID.get(chromPeak)].add(posPeak);
						numPeaks[condName2ID.get(condNames[i])]++;
					}
					i++;
				}//end of for loop
				
			}//end of else if condition
			
			else
				throw new IllegalArgumentException("Either you entered an invalid name" +
												   "for a format or this format is not supported yet.");
			
			// Assign peak lists to integer arrays and sort them 
			for(int i = 0; i < peakLists.length; i++)
				for(int j = 0; j < peakLists[i].length; j++) {
					peaks[i][j] = list2int(peakLists[i][j]);
					Arrays.sort(peaks[i][j]);
				}
			
			// Clear the lists as they are not useful anymore
			for(int i = 0; i < peakLists.length; i++)
				for(int j = 0; j < peakLists[i].length; j++)
					peakLists[i][j].clear();
			
			peakLists = null;
		}//end of try block
		
		catch(IOException e) {
			e.printStackTrace();
		}
		
		finally {
			try {
				for(File f:files) { 
					br = new BufferedReader(new FileReader(f)); 
					br.close(); 
				}
			} 
			catch (IOException e) { e.printStackTrace(); }
		}
	}//end of loadPeaks method
	
	
	public int[][][] getPeaks() { return peaks; }
	
	public Map<Integer,String> getID2hrom() { return id2chrom; }
	
	public Map<Integer, String> getID2condName() { return id2condName; }
	
	/**
	 * @return the number of peaks for each experiment (condition)
	 */
	public int[] getNumPeaks() { return numPeaks; }
	
	
	/**
	 * Takes as input a list of type <tt>Integer</tt> and returns an integer 
	 * array of primitive type
	 * @param list list of type <tt>Integer</tt>
	 * @return
	 */
	private int[] list2int(List<Integer> list) {
		int[] out = new int[list.size()];
		for(int i = 0; i < out.length; i++) {out[i] = list.get(i); }	   
		return out;
	}//end of list2int method


	/**>
	 * Usage: 					<br> 
	 * <pre>
	 * --file file1 --file file2 --file file3 ... --file fileN --format StatPeak --species "Mus musculus;mm8"
	 * </pre>
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		try {
			List<File> files = Args.parseFileHandles(args, "file");
			String format = Args.parseString(args, "format", "");
			Pair<Organism, Genome> organ_species = Args.parseGenome(args);
			Genome g = organ_species.cdr();
			
			if( (files.size() == 0) || format.equals("") || (organ_species == null) ) { printErrorMessage(); System.exit(-1); }
			
			PeakLoader pl = new PeakLoader(g, files, format);
			pl.loadPeaks();
			
			int foo = 3;
		} 
		
		catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}//end of main method
	
	
	private static void printErrorMessage() {
		System.err.println("Usage:\n " +
                "  --file <file in StatPeak or MACS format containing the peaks>\n" +
                "  (if it is more than one, just put multiple times the flag --file. E.g.: --file file1 --file file2 and so on)\n" +
                "  --format <MACS, StatPeak>\n" +
                "  (these are the two formats currently supported)\n" +
                "  --species <organism name;genome version>\n");
	}//end of printErrorMessage method

}//end of PeakLoader class
