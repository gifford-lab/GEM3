package edu.mit.csail.cgs.clustering.vectorcluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Vector;

/**
 * ElementParser creates a set of VectorClusterElement objects from a text
 * table in a file, with space separators.
 * 
 * @author Timothy Danford
 */
public class ElementParser {
	
	private File inputFile;
	private LinkedList<VectorClusterElement> elements;

	public ElementParser(File f) throws IOException {
		inputFile = f;
		parseFile();
	}
	
	private void parseFile() throws IOException {
		elements = new LinkedList<VectorClusterElement>();
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line;
		int lineCount = 0;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			String[] array = line.split("[\\s]+");
			Vector<Double> values = new Vector<Double>();
			for(int i = 0; i < array.length; i++) { 
				try { 
					double v = Double.parseDouble(array[i]);
					values.add(v);
				} catch(NumberFormatException nfe) { 
					values.add(null);
				}
			}
			DefaultVectorClusterElement de = new DefaultVectorClusterElement(values);
			de.addTag("line", String.valueOf(lineCount));
			lineCount+= 1;
			
			if(de.numMissingValues() < de.dimension()) {
				System.out.println("--> " + de.toString());
				elements.addLast(de); 
			}
		}
		br.close();
	}
	
	public Collection<VectorClusterElement> getElements() { return elements; }
}
