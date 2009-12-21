package edu.mit.csail.cgs.tools.motifs;

import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMetadataLoader;
import edu.mit.csail.cgs.datasets.chipchip.Experiment;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.BayesLocator;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.binding.BindingExpander;
import edu.mit.csail.cgs.ewok.verbs.BayesBindingGenerator;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.FastaWriter;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.ewok.verbs.RegionSorter;
import edu.mit.csail.cgs.utils.NotFoundException;

import java.util.*;
import java.io.*;
import java.sql.SQLException;

public class MotifFastaWriter {

	

	private static final double PROB_THRESHOLD = 0.3;
	private static final double SIZE_THRESHOLD = 2.0; //should be in range from 2.0 - 20 or 30

	private static final int PEAK_OFFSET_THRESHOLD = 200;
	private static final int PEAK_SEQUENCE_WINDOW_SIZE = 250;

	private static final String GENE_TABLE = "refGene";
	private static final int GENE_WINDOW_SIZE = 30000;	

	private static final String DELIM = ",";

	
	/**
	 * Written to read regions from a flat file with one region listed per line
	 * in the format:
	 * [chrom_number]	[start]-[end]
	 * e.g.
	 * 1	5043316-5047891
	 * 
	 * @param filename
	 * @return
	 */
	public static Vector<Region> getDataRegionsFromFile(Genome genome, String filename) {
		
		Vector<Region> dataRegions = new Vector<Region>();
		
		FileReader fileReader = null;
		BufferedReader reader = null;
		
		try {
			fileReader = new FileReader(filename);
			reader = new BufferedReader(fileReader);
			StreamTokenizer tokenizer = new StreamTokenizer(reader);			
			tokenizer.eolIsSignificant(true);
			
			int token = tokenizer.nextToken();
			while (token != StreamTokenizer.TT_EOF) {
				//read chromosome number
				String chrom;
				assert (token == StreamTokenizer.TT_NUMBER) || (token == StreamTokenizer.TT_WORD) : token;
				if (token == StreamTokenizer.TT_NUMBER) {
					int chromNumber = (int)tokenizer.nval;
					assert tokenizer.nval == (double)chromNumber;
					chrom = String.valueOf(chromNumber);
				}
				else {
					chrom = tokenizer.sval;					
				}
				
								
				//read start location
				token = tokenizer.nextToken();
				assert token == StreamTokenizer.TT_NUMBER;
				int start = (int)tokenizer.nval;
				assert tokenizer.nval == (double)start;
								
				/**
				 * read end location
				 * The hyphen will be interpreted as making the end location a 
				 * negative value, so check for that 
				 */
				token = tokenizer.nextToken();
				assert token == StreamTokenizer.TT_NUMBER;
				int end = (int)tokenizer.nval;
				assert tokenizer.nval == (double)end;
				assert end < 0;
				end = -end;
				
				//read EOL
				token = tokenizer.nextToken();
				assert (token == StreamTokenizer.TT_EOL) || (token == StreamTokenizer.TT_EOF);
				
				//System.out.println(tokenizer.lineno() + ": " + chrom + ", " + start + ", " + end);
				Region r = new Region(genome, chrom, start, end);
				dataRegions.add(r);
			
				//advance to next token
				token = tokenizer.nextToken();
			}
		}
		catch (IOException ioex) {
			ioex.printStackTrace();
		}
		finally {
			if (reader != null) {
				try {
					reader.close();
				}
				catch (IOException ioex2) {
					ioex2.printStackTrace();
				}
			}
		}
		
		return dataRegions;
	}
	
	
	/**
	 * 
	 * @param organism
	 * @param genome
	 * @param table
	 * @param expt
	 * @param version
	 * @param goodChroms
	 * @param filePrefix
	 * @param outputUnmatchedPeaks
	 * @param outputGenes
	 * @param outputNonGenePeaks
	 * @return
	 */
	public static Vector<Region> getDataRegionsFromPeaks(Organism organism, Genome genome, String table,
			String expt, String version, String type) {
		
		Vector<Region> dataRegions = new Vector<Region>();

		try {
			BindingScanLoader bsl = new BindingScanLoader();
			String bslVersion = expt + "," + version;
			Collection<BindingScan> bs = bsl.loadScans(genome, bslVersion, type);
			BindingExpander be = new BindingExpander(bsl, bs);
			
			//Scan through all the chromosomes
			Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				String chromName = currentChrom.getName();

				Vector<int[]> peaks = new Vector<int[]>();
				System.out.println("\n");
				//System.out.println("expt " + expts[i] + ", Chrom " + chromName + ", replicate " + j + ":");

				Iterator<BindingEvent> bindingIterator = be.execute(currentChrom);
				RegionSorter<BindingEvent> rs = new RegionSorter<BindingEvent>();
				Iterator<BindingEvent> sortedBindingIterator = rs.execute(bindingIterator); 
				if (!sortedBindingIterator.hasNext()) {
					System.out.println("expt " + expt + ", Chrom " + chromName + ": No Binding Events");
				}
				while (sortedBindingIterator.hasNext()) {
						BindingEvent event = sortedBindingIterator.next();
						int[] peak = new int[] {event.getStart(), event.getEnd()};
						peaks.addElement(peak);
				}
				System.out.println("expt " + expt + ", Chrom " + chromName + ": " + peaks.size() + " peaks");

				List<NamedRegion> peakRegions = processPeakData(peaks, genome, chromName, currentChrom.getWidth());					

				dataRegions.addAll(peakRegions);
			}
		}
		catch (SQLException sqlex) {
			sqlex.printStackTrace();
		}


		return dataRegions;
	}
	
	
	/**
	 * Combine peaks that are closest to each other and generate regions
	 * @param repPeakInfo
	 * @return
	 */
	public static List<NamedRegion> processPeakData(Vector<int[]> peaks, Genome genome, String chromName, int chromLength) {
		List<NamedRegion> processedPeaks = new Vector<NamedRegion>();
		int index = 0;

		Vector<int[]> peakSequence = findNearbyPeakSequence(peaks, index, PEAK_OFFSET_THRESHOLD);
		while (peakSequence.size() > 0) {
			//System.out.println("Sequential peaks: " + peakSequence.size());
			int start = (peakSequence.elementAt(0))[0];
			int end = (peakSequence.elementAt(peakSequence.size()-1))[1];
			
			String regionName = chromName + ":" + start + "-" + end;
			
			int middle = (start + end) / 2;
			Gene[] closestGenes = findClosestGenes(genome, GENE_TABLE, chromName, middle, GENE_WINDOW_SIZE, GENE_WINDOW_SIZE);
			
			if ((closestGenes[0] != null) || (closestGenes[1] != null)) {
				int upDist = Integer.MAX_VALUE;
				int downDist = Integer.MAX_VALUE;
				if (closestGenes[0] != null) {
					upDist = start - closestGenes[0].getStart();  
				}
				if (closestGenes[1] != null) {
					downDist = closestGenes[1].getStart() - start; 
				}
			  
				if ( upDist < downDist) {
					regionName = regionName + ", " + downDist + " bp downstream of " + closestGenes[0].getName(); 
				}
				else {
					regionName = regionName + ", " + upDist + " bp upstream of " + closestGenes[1].getName();
				}
			}
			int sequenceStart = (int)Math.max(start - (PEAK_SEQUENCE_WINDOW_SIZE/2), 0);
			int sequenceEnd = (int)Math.min(end + (PEAK_SEQUENCE_WINDOW_SIZE/2), chromLength);
			NamedRegion region = new NamedRegion(genome, chromName, sequenceStart, sequenceEnd, regionName);
			processedPeaks.add(region);
			
			index = index + peakSequence.size();
			peakSequence = findNearbyPeakSequence(peaks, index, PEAK_OFFSET_THRESHOLD);			
		}

		return processedPeaks;
	}
	
	

	
	
	/**
	 * 
	 * @param dataRegions
	 * @param outputFilename
	 */
	public static void writeDataRegions(Vector<Region> dataRegions, String outputFilename) {
		
		FastaWriter writer = null;
		try {
			writer = new FastaWriter(outputFilename);
			writer.consume(dataRegions.iterator());
		}
		catch (FileNotFoundException fnfex) {
			fnfex.printStackTrace();
		}
		catch (IOException ioex) {
			ioex.printStackTrace();			
		}
		finally {
			if (writer != null) {
				writer.close();
			}
		}
	}
	
	
	/**
	 * Find a sequence of consecutive nearby peaks
	 * This method should eventually move into another class
	 *  
	 * @param repPeakInfo
	 * @param indices
	 * @param maxDistance
	 * @return
	 */
	public static Vector<int[]> findNearbyPeakSequence(Vector<int[]> peaks, int index, int maxDistance) {
		Vector<int[]> peakSequence = new Vector<int[]>();

		int prevPeakLocation = -1;

		//find the first peak from among the replicates
		if (index < peaks.size()) {
			prevPeakLocation = (peaks.elementAt(index))[0];
		}
		
		if (prevPeakLocation == -1) {
			//no peaks, return an empty sequence
			return peakSequence;
		}
		else {
			peakSequence.addElement(peaks.elementAt(index));
			index++;	
		}


		/**
		 * continue finding peaks until consecutive peaks are farther apart than
		 * the specified max distance 
		 */		
		boolean done = false;
		while (!done) {
			int nextPeakLocation = -1;
			if (index < peaks.size()) {		
				nextPeakLocation = (peaks.elementAt(index))[0];
			}
			if ((nextPeakLocation == -1) || ((nextPeakLocation - prevPeakLocation) > maxDistance)) {
				/**
				 * stop if no more peaks could be found, or consecutive peaks
				 * are too far apart 
				 */
				done = true;
			}
			else { 
				//peaks are close together
				peakSequence.addElement(peaks.elementAt(index));
				index++;
				prevPeakLocation = nextPeakLocation;
			}
		}

		return peakSequence;
	}

	
	/**
	 * Find the genes closest to the specified location
	 * This method should eventually move into another class
	 * @param genome
	 * @param chromName
	 * @param loc
	 * @return
	 */
	public static Gene[] findClosestGenes(Genome genome, String table, String chromName, int loc, int upLimit, int downLimit) {
		RefGeneGenerator rgg = new RefGeneGenerator(genome, table);

		int regionStart = (int)Math.max(loc - upLimit, 0);
		int regionEnd = (int)(loc + downLimit);
		Region region = new Region(genome, chromName, regionStart, regionEnd);
		Iterator<Gene> geneIter = rgg.execute(region);

		Gene[] closestGenes = new Gene[2];
		Arrays.fill(closestGenes, null);
		int[] closestDistances = new int[] {Integer.MAX_VALUE, Integer.MAX_VALUE };
		while (geneIter.hasNext()) {
			Gene currentGene = geneIter.next();
			int currentStart = currentGene.getStart();
			int currentDist = (int)Math.abs(loc - currentStart);
			if (currentStart < loc) {				
				if (currentDist < closestDistances[0]) {
					closestGenes[0] = currentGene;
					closestDistances[0] = currentDist;
				}	
			}
			else {
				if (currentDist < closestDistances[1]) {
					closestGenes[1] = currentGene;
					closestDistances[1] = currentDist;					
				}
			}			
		}

		return closestGenes;
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/**
		 * write out fasta files for the regions of the genome tiled by the
		 * Hox arrays in the ppg project
		 */
		
//		try {
//			String filename = "G:\\projects\\ppg\\domains\\well_tiled_regions.txt";
//			String outputFilename = "";
//
//			Organism mouse = Organism.getOrganism("Mus musculus");
//			Genome mm8 = mouse.getGenome("mm8");
//			Vector<Region> dataRegions = MotifFastaWriter.getDataRegionsFromFile(mm8, filename);
//			MotifFastaWriter.writeDataRegions(dataRegions, outputFilename);
//		} 
//		catch (NotFoundException ex) {
//			// TODO Auto-generated catch block
//			ex.printStackTrace();
//		}
		
		/**
		 * write out fasta files for regions where peaks have been called
		 */
		try {
			Organism mouse = Organism.getOrganism("Mus musculus");
			Genome mm8 = mouse.getGenome("mm8");
			String expt = "Mm Hb9:HBG3:Hb9 Stage vs WCE:HBG3:Hb9 Stage";
			String version = "1/25/07, default params";
			String outputFilename = "hb9_peak_sequences.fasta";
			String type = "BayesBindingGenerator";
			Vector<Region> dataRegions = MotifFastaWriter.getDataRegionsFromPeaks(mouse, mm8, GENE_TABLE, expt, version, type);
			MotifFastaWriter.writeDataRegions(dataRegions, outputFilename);
		} 
		catch (NotFoundException ex) {
			// TODO Auto-generated catch block
			ex.printStackTrace();
		}
	}

}
