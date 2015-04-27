package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqAlignment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.ExtReadHit;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Loads reads from a file. Formats supported:
 * ELAND, NOVO, BOWTIE, BED								<br>
 * 
 * <tt>NOTE</tt>: For each new file loader when we create the method
 * <tt>countReads()</tt>,											<br> 
 * we should always call at the end of this method
 * the <tt>AlignmentFileReader.populateArrays</tt>, 				<br>
 * since it not only populates the desired entries but sorts them out as well, 
 * something that is necessary.
 * @author shaun
 *
 */
public class FileReadLoader extends ReadLoader{

	protected List<File> files=null;
	protected String format;
	// a list of file readers, each for one file (replicate)
	protected List<AlignmentFileReader> fileReaders = new ArrayList<AlignmentFileReader>();
	protected List<PairedAlignmentFileReader> pairedReaders = new ArrayList<PairedAlignmentFileReader>();
	protected int maxMismatch=0;
	protected boolean useNonUnique=true;
	protected int currID=0;

	public FileReadLoader(List<File> f, String format){
		// estimate genome chrom length by combining all files
		for(File file : f){
			if(!file.isFile()){System.err.println("File not found: "+file.getName());System.exit(1);}
			if(format.equals("SAM")){
				SAMReader currReader = new SAMReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("TOPSAM")){
				TophatSAMReader currReader = new TophatSAMReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("ELAND")){
				ElandFileReader currReader = new ElandFileReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("NOVO")){
				NovoFileReader currReader = new NovoFileReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("BOWTIE")){
				BowtieFileReader currReader = new BowtieFileReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("BED")){
				BEDFileReader currReader = new BEDFileReader(file);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("BEDPE")){
                BEDPEFileReader currReader = new BEDPEFileReader(file);
                pairedReaders.add(currReader);
                currID = currReader.getCurrID();
			}else{
			    System.err.println("Unknown file format: "+format);
			    System.exit(1);
			}
		}

		//Combine the chromosome information
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		for(AlignmentFileReader a : fileReaders){
			Map<String, Integer> currMap = a.getGenome().getChromLengthMap();
			for(String s: currMap.keySet()){
				if(!chrLenMap.containsKey(s) || chrLenMap.get(s)<currMap.get(s))
					chrLenMap.put(s, currMap.get(s));
			}
		}
		gen = new Genome("Genome", chrLenMap);
	}

	public FileReadLoader(Genome g, List<File> f, String format, int maxMismatch, boolean useNonUnique, int rLen, int idSeed){
		super(g, rLen);
		
		this.maxMismatch=maxMismatch;
		this.useNonUnique=useNonUnique;
		files=f;
		if(format==null){this.format="BED";}
		else{this.format=format;}
		currID = idSeed;
			
		for(File file : files){
			if(!file.isFile()){System.err.println("File not found: "+file.getName());System.exit(1);}
			if(format.equals("SAM")){
				SAMReader currReader = new SAMReader(file,gen,maxMismatch,useNonUnique, currID, chrom2ID, id2Chrom);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("TOPSAM")){
				TophatSAMReader currReader = new TophatSAMReader(file,gen,maxMismatch,useNonUnique, currID, chrom2ID, id2Chrom);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("ELAND")){
				ElandFileReader currReader = new ElandFileReader(file,gen,maxMismatch,useNonUnique, currID, chrom2ID, id2Chrom);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("NOVO")){
				NovoFileReader currReader = new NovoFileReader(file,gen,useNonUnique, currID, chrom2ID, id2Chrom);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("BOWTIE")){
				BowtieFileReader currReader = new BowtieFileReader(file,gen,useNonUnique, currID, chrom2ID, id2Chrom);
				fileReaders.add(currReader);
				currID = currReader.getCurrID();
			}else if(format.equals("BED")){
                BEDFileReader currReader = new BEDFileReader(file,gen,useNonUnique, currID, chrom2ID, id2Chrom);
                fileReaders.add(currReader);
                currID = currReader.getCurrID();
            }else if(format.equals("BEDPE")){
                BEDPEFileReader currReader = new BEDPEFileReader(file,gen, -1, useNonUnique, currID, chrom2ID, id2Chrom);
                pairedReaders.add(currReader);
                currID = currReader.getCurrID();
            }else{
			    System.err.println("Unknown file format: "+format);
			    System.exit(1);
			}
		}

		for(AlignmentFileReader a : fileReaders){
			totalHits+=a.getTotalHits();
			totalWeight+=a.getTotalWeight();
			//Set the read length to the maximum value for now, fix later 
			if(a.getReadLength()>readLength)
				readLength = a.getReadLength();
		}		
	}

	public int getCurrID(){return currID;}
	
	//Since we counted the reads as we added them, we don't need to explicitly count again
	protected double countHits() {
		return(totalHits);
	}

	//Load the reads from our files
	public List<ReadHit> loadHits(Region r) {
		ArrayList<ReadHit> hits = new ArrayList<ReadHit>();
		for(AlignmentFileReader a : fileReaders){
			hits.addAll(a.loadHits(r));
		}
		return hits;
	}
	
	//Load paired reads from our files
	public List<ReadHit> loadPairs(Region r) {
		ArrayList<ReadHit> hits = new ArrayList<ReadHit>();
		//TODO: NOT YET IMPLEMENTED FOR FILES
		return hits;
	}
	//Load paired reads from our files
	public List<PairedHit> loadPairsAsPairs(Region r) {
		ArrayList<PairedHit> hits = new ArrayList<PairedHit>();
		//TODO: NOT YET IMPLEMENTED FOR FILES
		return hits;
	}
	
	// load all start coordinates (unsorted if multiple files)
	public ArrayList<int [][][]> getAllStarts(){
		ArrayList<int [][][]> allStarts = new ArrayList<int [][][]>();
		for(AlignmentFileReader a : fileReaders){
			allStarts.add(a.getFivePrimes());
		}
		return allStarts;
	}
	
	// load paired read hit coordinates (sorted) and counts
	public Pair<ArrayList<Integer>,ArrayList<Float>> loadStrandedFivePrimeCounts(Region r, char strand){
		//  TreeMap<sorted key: 5' coordinates, value: counts>
		// The duplicate hits will be consolidated as bases and counts using TreeMap
		TreeMap<Integer,Float> allHits = new TreeMap<Integer,Float>();
		for(AlignmentFileReader a : fileReaders){
			Pair<ArrayList<Integer>,ArrayList<Float>> bases = a.loadStrandedFivePrimeCounts(r, strand);
			ArrayList<Integer> coordinates = bases.car();
			ArrayList<Float> counts = bases.cdr();
			for (int i=0;i<coordinates.size();i++){
				int coordinate = coordinates.get(i);
	        	if (allHits.containsKey(coordinate))
	        		allHits.put(coordinate, allHits.get(coordinate)+counts.get(i));
	        	else
	        		allHits.put(coordinate, counts.get(i));	// hit with same coordinate should have same weight
			}
		}
		ArrayList<Integer> coords = new ArrayList<Integer>();
		coords.addAll(allHits.keySet());
		ArrayList<Float> counts = new ArrayList<Float>();
		counts.addAll(allHits.values());
		return new Pair<ArrayList<Integer>,ArrayList<Float>>(coords, counts);
	}

	// load coordinates and counts for paired-end data
	// cheats for now, assumes single file only
	public Pair<Pair<ArrayList<Integer>, ArrayList<ArrayList<Integer>>>,ArrayList<ArrayList<Float>>> loadStrandedPairedEndCounts(Region r, char strand){
	    PairedAlignmentFileReader a = pairedReaders.get(0);
	    return a.loadStrandedFivePrimeCounts(r, strand);

	}
	
	//Load the extended reads from our files
	public List<ExtReadHit> loadExtHits(Region r, int startShift, int fivePrimeExt, int threePrimeExt) {
		ArrayList<ExtReadHit> hits = new ArrayList<ExtReadHit>();
		for(AlignmentFileReader a : fileReaders){
			List<ReadHit> tmp = a.loadHits(r);
			for(ReadHit h : tmp)
				hits.add(new ExtReadHit(h,startShift, fivePrimeExt, threePrimeExt));
			/*List<Read> reads = a.loadReads(r);
			//Filter
			for(Read x : reads){
				for(ReadHit h : x.getHits()){
					if(h.overlaps(r))
						hits.add(new ExtReadHit(h,startShift, fivePrimeExt, threePrimeExt));					
				}
			}*/
		}
		return hits;
	}
	
	// count number of reads in region
	public int countHits(Region r) {
		int count = 0;
		for(AlignmentFileReader a : fileReaders){
			count += a.loadHits(r).size();
		}
		return count;
	}
	// sum weights of reads in region
	public double sumWeights(Region r) {
		double sum = 0;
		for(AlignmentFileReader a : fileReaders){
			for(ReadHit h : a.loadHits(r))
				sum += h.getWeight();
		}
		return sum;
	}
	//Stranded read hit count
	protected double countStrandedWeight(char strand) {
		double count=0;
		for(AlignmentFileReader a : fileReaders){
			count+=a.getStrandedTotalWeight(strand);
		}
		return count;
	}
	//Set a new genome
	public void setGenome(Genome g){
		gen=g;
		for(AlignmentFileReader a : fileReaders){
			a.setGenome(g);
		}
	}
	
	public void cleanup(){
		for(AlignmentFileReader a : fileReaders){
			a.cleanup();
		}
	}
}

