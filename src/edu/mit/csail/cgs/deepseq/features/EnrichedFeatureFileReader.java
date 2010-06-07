package edu.mit.csail.cgs.deepseq.features;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Genome;

public class EnrichedFeatureFileReader {
	Genome gen;
	ArrayList<EnrichedFeature> features = new ArrayList<EnrichedFeature>();
	
	public EnrichedFeatureFileReader(Genome g){}
	public EnrichedFeatureFileReader(Genome g, String fileName){
		execute(g, fileName);
	}
	
	//Accessor
	public ArrayList<EnrichedFeature> getFeatures(){return features;}
	
	//Read a file
	public ArrayList<EnrichedFeature> execute(Genome g, String fileName){
		File f = new File(fileName);
		return(execute(g, f));		
	}
	public ArrayList<EnrichedFeature> execute(Genome g, File f){
		features = new ArrayList<EnrichedFeature>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line = reader.readLine(); //remove header
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	            String[] words = line.split("\\t");
	            String[] coordsA = words[0].split(":");
	            String[] coordsB = coordsA[1].split("-");
	            String cchr=coordsA[0]; char strand = coordsA.length>=3 ? coordsA[2].charAt(0) : '.';
	            int cstart=new Integer(coordsB[0]), cend=new Integer(coordsB[1]);
	            String[] peak = words[2].split(":");
	            String pchr=peak[0]; 
	            int pstart=new Integer(peak[1]);
	            
	            EnrichedFeature e = new EnrichedFeature();
	            e.coords=new Region(g, cchr, cstart, cend);
	            e.peak=new Point(g, pchr, pstart);
	            e.signalMaxHits = new Double(words[4]);
	            e.backMaxHits = new Double(words[5]);
	            e.signalTotalHits = new Double(words[7]);
	            e.backTotalHits = new Double(words[8]);
	            e.score = new Double(words[6]);
	            e.overrep = new Double(words[9]);
	            e.strand = strand;
	            	
	            features.add(e);
			}reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(features);
	}
}
