/*
 * Author: tdanford
 * Date: Aug 24, 2008
 */
package edu.mit.csail.cgs.datasets.alignments.parsing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parses a CLUSTALW output file into separate gapped alignments.  
 * 
 * Takes an optional regexp parameter, as a way to pull out an optional tag  
 * from the clustal sequence names.  
 * 
 * Indexes the internal strings by the species tag (or by the clustal sequence name, 
 * if no tag regexp is given).  
 * 
 * @author tdanford
 *
 */
public class ClustalParser { 
	
	/**
	 * A pattern used for recognizing the tags/names of sequences</br>
	 * Default value is: [\\S| ]+</br>
	 * Namely, it allows any non-white space character plus an one-spaced blank
	 */
	private Pattern clustalSpeciesPatt; 
	private String lineBreakPatt;
	
	/**
	 * Map which contains as keys the tags/names of the sequences and 
	 * as their corresponding values the sequences themselves
	 */
	private Map<String,StringBuilder> sequences;
	
	public ClustalParser() { 
		sequences = new HashMap<String,StringBuilder>();
		clustalSpeciesPatt = Pattern.compile("(.*)");
		lineBreakPatt = "\\t|\\s{2,}";
	}
	
	public ClustalParser(Pattern speciesPattern) { 
		this();
		clustalSpeciesPatt = speciesPattern;
	}
	
	public ClustalParser(Pattern speciesPattern, String lineBreakPattern) { 
		this();
		clustalSpeciesPatt = speciesPattern;
		lineBreakPatt = lineBreakPattern;
	}
	
	public ClustalParser(File f) throws IOException { 
		this();
		parseFile(f);
	}
	
	public ClustalParser(File f, Pattern speciesPattern) throws IOException { 
		this(speciesPattern);
		parseFile(f);
	}
	
	public ClustalParser(File f, Pattern speciesPattern, String linePatt) throws IOException { 
		this(speciesPattern, linePatt);
		parseFile(f);
	}
	
	public String getSequence(String species) { return sequences.get(species).toString(); }
	public Set<String> getSpecies() { return sequences.keySet(); }
	
	public void parseFile(File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null) { 
			parseLine(line);
		}
		br.close();
	}
	
	/**
	 * Scans a file for species names ('tags'), without actually parsing the file's 
	 * sequence itself.  
	 * 
	 * @param f
	 * @return
	 * @throws IOException
	 */
	public Set<String> speciesScan(File f) throws IOException {
		TreeSet<String> species = new TreeSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null) { 
			String[] array = line.trim().split(lineBreakPatt);
			if(array.length==2) {
				String tag = array[0];
				if(clustalSpeciesPatt != null) { 
					Matcher m = clustalSpeciesPatt.matcher(array[0]);
					if(m.matches()) { 
						tag = m.group(1); 
					} else { 
						tag = null; 
					}
				}				

				if(tag != null) { 
					species.add(tag); 
				}
			}
		}
		br.close();		
		return species;
	}
	
	/**
	 * It parses the line which supposedly contains the species name as well as the sequence
	 * itself and possibly some number indexing.</br>
	 * Here we make the strong assumption that between the species name and the sequences there is 
	 * at least a two-blank interval or a tab delimiter and that the sequence species names can't have
	 * more than one-spaced blank in between.
	 * @param line
	 * @throws IOException
	 */
	public void parseLine(String line) throws IOException
	{ 
		
		if( !line.matches("[*.:\\s]+") )
		{
			String[] array = line.trim().split(lineBreakPatt);
			
			// There is no indexing at the end of each sequence
			if(array.length == 2)
			{
				String[] purifySeq = array[1].trim().split("\\s+");
				array[1] = purifySeq[0];
				populateSequencesMap(array, clustalSpeciesPatt, sequences);
			}
			
			// There is indexing at the end of each sequence
			else if(array.length == 3)
			{
				// If the last token is not a number indicating the current length of the sequence
				// read so far, an exception will be thrown.
				try 
				{
					int index = Integer.parseInt(array[2]);
				} 
				catch(NumberFormatException nfe) 
				{
					throw new IOException(String.format("Couldn't parse element: %s", array[2]));
				}
				
				String[] purifySeq = array[1].trim().split("\\s+");
				array[1] = purifySeq[0];
				populateSequencesMap(array, clustalSpeciesPatt, sequences);
			}
		}
	
	}//end of parseLine method
	
	private void populateSequencesMap(String[] array, Pattern clustalSpeciesPatt, Map<String,StringBuilder> seqMap)
	{
		String tag = array[0];
		
		//System.out.println(String.format("PSM: %s \"%s\"", array[0], array[1]));
		
		if(clustalSpeciesPatt != null) { 
			Matcher m = clustalSpeciesPatt.matcher(array[0]);
			if(m.matches()) {
				tag = m.group(1);
			} else { 
				tag = null; 
			}
		}
		
		if(tag != null) { 
			String seq = array[1];
			if(!seqMap.containsKey(tag)) { 
				seqMap.put(tag, new StringBuilder());
			}
			seqMap.get(tag).append(seq);
		}
	}//end of populateSequencesMap method
}
