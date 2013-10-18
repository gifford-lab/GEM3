/**
 * 
 */
package edu.mit.csail.cgs.deepseq.utilities;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.imageio.ImageIO;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixPainter;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.strings.StringUtils;

/**
 * @author Yuchun Guo
 * This clsss is to hold some frequently used static methods
 *
 */
public class CommonUtils {
	/*
	 * load text file in CGS Point format
	 * chr:coord, e.g. 1:234234
	 */
	public static ArrayList<Point> loadCgsPointFile(String filename, Genome genome) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				String f[] = line.trim().split("\t");
				if (f[0].startsWith("#"))
					continue;
				if (f[0].startsWith("Position"))		// to be compatible with GPS/GEM event file
					continue;
				Region point = Region.fromString(genome, f[0]);
				if (point!=null)
					points.add(new Point(genome, point.getChrom(),point.getStart()));
			}
		}
		catch(IOException ioex) {
			ioex.printStackTrace();
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
		return points;
	}
	
	public static ArrayList<Region> loadRegionFile(String fname, Genome gen){
		ArrayList<Region> rset = new ArrayList<Region>();
		try{
			File rFile = new File(fname);
			if(!rFile.isFile()){
				System.err.println("\nThe region file is not found: "+fname+"!");
				System.exit(1);
			}
	        BufferedReader reader = new BufferedReader(new FileReader(rFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            if (line.startsWith("#"))
					continue;
	            String[] words = line.split("\\s+");
            	Region r = Region.fromString(gen, words[0]);
            	if (r!=null)
            		rset.add(r);
            }
	        reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return rset;
	}
	
	/**
	 * Load ENCODE narrow Peak format, sorted by score, then by signal<br>
	 * 
	 *  Col1: chromosome name 
		Col2: start position of peak region (0-based)
		Col3: stop position of peak region (1-based)
		Col4: peak name (any string)
		Col5: peak score (any scoring/ranking measure for the peaks)
		Col6: Strand (Always a '.' i.e. no orientation)
		Col7: signalValue (any scoring/ranking measure for the peaks)
		Col8: -log10(pValue) if the peak caller does not provide a pValue set this Column to -1
		Col9: -log10(qValue) if	the peak caller	does not provide a qValue set this Column to -1
		Col10: peak summit (0-based offset from start i.e. col2) single base best prediction of binding event/peak summit
	 * @param genome
	 * @param filename
	 * @return
	 */
	static public ArrayList<NarrowPeak> load_narrowPeak(Genome genome, String filename, boolean isSorted) {
		CommonUtils util = new CommonUtils();
		ArrayList<NarrowPeak> results = new ArrayList<NarrowPeak>();
		ArrayList<String> txt = readTextFile(filename);
		for (String s:txt){
			if (s.startsWith("#"))
				continue;
			String[] f = s.split("\t");
			for (int i=0;i<f.length;i++){
				if (f[i].equalsIgnoreCase("Inf"))
					f[i]="999";
			}
			NarrowPeak p = util.new NarrowPeak(genome, f[0], Integer.parseInt(f[1]), Integer.parseInt(f[2]), Double.parseDouble(f[4]), Double.parseDouble(f[6]), Double.parseDouble(f[7]), Double.parseDouble(f[8]), Integer.parseInt(f[9]));
			results.add(p);
		}
		results.trimToSize();
		if (!isSorted)
			Collections.sort(results);
		return results;
	}
	
	public class NarrowPeak implements Comparable<NarrowPeak> {
		public Region region;
		public Point summit;
		public double score;
		public double signal;
		public double pvalue;
		public double qvalue;
		public NarrowPeak(Genome genome, String chrom, int start, int end, double score, double signal, double pvalue, double qvalue, int summitOffset){
			String chr = chrom.replace("chr", "").replace("Chr", "");
			region = new Region(genome, chr, start, end);
			summit = new Point(genome, chr, start+summitOffset);
			this.score = score;
			this.signal = signal;
			this.pvalue = pvalue;
			this.qvalue = qvalue;
		}
		@Override
		public int compareTo(NarrowPeak f) {
			if(score>f.score){return(-1);}
			else if(score<f.score){return(1);}
			else return(0);
		}
		public int compareToBySignal(NarrowPeak f) {
			if(signal>f.signal){return(-1);}
			else if(signal<f.signal){return(1);}
			else return(0);
		}
	}
	
	/** load text file in SISSRS output BED format, then sort by p-value
	 * 	Chr		cStart	cEnd	NumTags	Fold	p-value
		---		------	----	-------	----	-------
		chr1	4132791	4132851	38		71.44	5.0e-006
	 */
	static public ArrayList<SISSRS_Event> load_SISSRs_events(Genome genome, String filename, boolean isPreSorted) {

		CommonUtils util = new CommonUtils();
		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<SISSRS_Event> events = new ArrayList<SISSRS_Event>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			// skip header, data starts at 58th line
			for (int i=1;i<58;i++){
				bin.readLine();
			}
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				String[] t = line.split("\\t");
				if (t.length==6){	// region format
					events.add(util.new SISSRS_Event(genome, t[0].replaceFirst("chr", ""), Integer.parseInt(t[1]), Integer.parseInt(t[2]),
						Double.parseDouble(t[3]), Double.parseDouble(t[4]), Double.parseDouble(t[5])));
				}
				if (t.length==5){	// point format
					events.add(util.new SISSRS_Event(genome, t[0].replaceFirst("chr", ""), Integer.parseInt(t[1]), Integer.parseInt(t[1]),
						Double.parseDouble(t[2]), Double.parseDouble(t[3]), Double.parseDouble(t[4])));
				}
			}
		}
		catch(IOException ioex) {
			ioex.printStackTrace();
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
		// sort by pvalue
		if (!isPreSorted)
			Collections.sort(events);
		
		return events;
	}
	
	public class SISSRS_Event implements Comparable<SISSRS_Event>{
		public Region region;
		public double tags;
		public double fold;
		public double pvalue;
		public SISSRS_Event(Genome g, String chrom, int start, int end, double tags, double fold, double pvalue){
			this.region = new Region(g, chrom, start, end);
			this.tags = tags;
			this.fold = fold;
			this.pvalue = pvalue;
		}
		public Point getPeak(){
			return region.getMidpoint();
		}
		//Comparable default method
		public int compareTo(SISSRS_Event p) {
			double diff = pvalue-p.pvalue;
			return diff==0?0:(diff<0)?-1:1;
		}
	}
	
	public static String timeElapsed(long tic){
		return timeString(System.currentTimeMillis()-tic);
	}
	private static String timeString(long length){
		float sec = length/1000F;
		return sec>3600?
				String.format("%.1fh",sec/3600):
				sec>60?
					String.format("%.1fm",sec/60):
					sec>1?
						String.format("%.1fs",sec):
						String.format("%dmils",length)	;
	}
	public static String getDateTime() {
        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
        Date date = new Date();
        return dateFormat.format(date);
    }
	public static void printArray(double[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static void printArray(int[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static String arrayToString(int[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%d\t",array[i]));
        }
        output.append(String.format("%d",array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%.2f\t",array[i]));
        }
        output.append(String.format("%.2f",array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array, String format){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format(format+"\t",array[i]));
        }
        output.append(String.format(format,array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array, int digit){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%."+digit+"f\t",array[i]));
        }
        output.append(String.format("%."+digit+"f",array[array.length-1]));
        return output.toString();
	}	
	public static String matrixToString(double[][] matrix, int digit, String[] names){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<matrix.length;i++){
        	if (names!=null)
        		output.append(names[i]).append("\t");
        	output.append(arrayToString(matrix[i], digit)).append("\n");
        }
        return output.toString();
	}
	public static void replaceEnd(StringBuilder sb, char newChar){
		sb.deleteCharAt(sb.length()-1).append(newChar);
	}
	
	public static void writeFile(String fileName, String text){
		try{
			FileWriter fw = new FileWriter(fileName, false); //new file
			fw.write(text);
			fw.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
//		System.out.println("File was written to "+fileName);
	}
	
	public static ArrayList<String> readTextFile(String fileName){
		ArrayList<String> strs = new ArrayList<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
	        String line;
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            if (line.length()!=0)
	            	strs.add(line);
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+fileName);
            e.printStackTrace(System.err);
        }   
        return strs;
	}
	
	public static ArrayList<String> readFastaFile(String fileName){
		ArrayList<String> strs = new ArrayList<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
	        String line;
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            if (!line.startsWith(">") && line.length()!=0)
	            	strs.add(line);
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+fileName);
            e.printStackTrace(System.err);
        }   
        return strs;
	}
	
	/** 
	 * Find the index that gives value larger than or equal to the key
	 * @param values
	 * @param key
	 * @return
	 */
	public static int findKey(double[] values, double key){
		int index = Arrays.binarySearch(values, key);
		if( index < 0 )  							// if key not found
			index = -(index+1); 
		else{		// if key match found, continue to search ( binarySearch() give undefined index with multiple matches)
			for (int k=index; k>=0;k--){
				if (values[k]==key)
					index = k;
				else
					break;
			}
		}
		return index;
	}
	
	public static String padding(int repeatNum, char padChar) throws IndexOutOfBoundsException {
	      if (repeatNum < 0) {
	          throw new IndexOutOfBoundsException("Cannot pad a negative amount: " + repeatNum);
	      }
	      if (repeatNum == 0) 
	    	  return "";

	      final char[] buf = new char[repeatNum];
	      for (int i = 0; i < buf.length; i++) {
	          buf[i] = padChar;
	      }
	      return new String(buf);
	}
	public static String padding(int repeatNum, String s) throws IndexOutOfBoundsException {
	      if (repeatNum < 0) {
	          throw new IndexOutOfBoundsException("Cannot pad a negative amount: " + repeatNum);
	      }
	      if (repeatNum == 0) 
	    	  return "";
	      
	      StringBuilder sb = new StringBuilder();
	      for (int i = 0; i < repeatNum; i++) {
	          sb.append(s);
	      }
	      return sb.toString();
	}
    /*
     * parse motif string, version, threshold, species_id from command line
     */
    public static Pair<WeightMatrix, Double> loadPWM(String args[], int orgId){
    	Pair<WeightMatrix, Double> pair=null;
    	try {   
    		double motifThreshold = Args.parseDouble(args, "motifThreshold", -1);
    		if (motifThreshold==-1){
    			System.err.println("No motif threshold was provided, default=9.99 is used.");
    			motifThreshold = 9.99;
    		}  
      	    String motifString = Args.parseString(args, "motif", null);
            if (motifString!=null){
      	      String motifVersion = Args.parseString(args, "version", null);
      	      
      		  int motif_species_id = Args.parseInteger(args, "motif_species_id", -1);
//      		  int wmid = WeightMatrix.getWeightMatrixID(motif_species_id!=-1?motif_species_id:orgId, motifString, motifVersion);
      		  int wmid = WeightMatrix.getWeightMatrixID(motifString, motifVersion);
      		  WeightMatrix motif = WeightMatrix.getWeightMatrix(wmid);
      		  System.out.println(motif.toString());
      		  pair = new Pair<WeightMatrix, Double>(motif, motifThreshold);
            }
            else{
            	String pfmFile = Args.parseString(args, "pfm", null);
            	if (pfmFile!=null){
            		double gc = Args.parseDouble(args, "gc", 0.41);
            		WeightMatrix motif = loadPWM_PFM_file(pfmFile, gc);
                	pair = new Pair<WeightMatrix, Double>(motif, motifThreshold);
            	}
            }
          } 
          catch (NotFoundException e) {
            e.printStackTrace();
          }  
  		  return pair;
    }
    /**
     * Load single (first) PWM from PFM files
     * @param pfmFile	PFM file in STAMP format (simplified TRANSFAC format)
     * @param gc	expected gc fraction
     * @return
     */
    public static WeightMatrix loadPWM_PFM_file(String pfmFile, double gc){
		try{
			List<WeightMatrix> wms = WeightMatrixImport.readTRANSFACFreqMatrices(pfmFile, "file");
			if (wms.isEmpty()){
				return null;
			}
			else{		// if we have valid PFM
				WeightMatrix wm = wms.get(0);
				float[][] matrix = wm.matrix;
				// normalize
		        for (int position = 0; position < matrix.length; position++) {
		            double sum = 0;
		            for (int j = 0; j < WeightMatrix.letters.length; j++) {
		                sum += matrix[position][WeightMatrix.letters[j]];
		            }
		            for (int j = 0; j < WeightMatrix.letters.length; j++) {
		                matrix[position][WeightMatrix.letters[j]] = (float)(matrix[position][WeightMatrix.letters[j]] / sum);
		            }
		        }
		        // log-odds
		        for (int pos = 0; pos < matrix.length; pos++) {
		            for (int j = 0; j < WeightMatrix.letters.length; j++) {
		                matrix[pos][WeightMatrix.letters[j]] = (float)Math.log(Math.max(matrix[pos][WeightMatrix.letters[j]], .001) / 
		                		(WeightMatrix.letters[j]=='G'||WeightMatrix.letters[j]=='C'?gc/2:(1-gc)/2));
		            }
		        } 
		        return wm;
			}
		}
		catch (IOException e){
			return null;
		}
    }
 
    /**
     * Load multiple PWMs from PFM files
     * @param pfmFile	PFM file in STAMP format (simplified TRANSFAC format)
     * @param gc	expected gc fraction
     * @return
     */
    public static List<WeightMatrix> loadPWMs_PFM_file(String pfmFile, double gc){
		try{
			List<WeightMatrix> wms = WeightMatrixImport.readTRANSFACFreqMatrices(pfmFile, "file");
			if (wms.isEmpty()){
				return null;
			}
			else{		// if we have valid PFM
				for (int i=0;i<wms.size();i++){
					WeightMatrix wm = wms.get(i);
					float[][] matrix = wm.matrix;
					// normalize
			        for (int position = 0; position < matrix.length; position++) {
			            double sum = 0;
			            for (int j = 0; j < WeightMatrix.letters.length; j++) {
			                sum += matrix[position][WeightMatrix.letters[j]];
			            }
			            for (int j = 0; j < WeightMatrix.letters.length; j++) {
			                matrix[position][WeightMatrix.letters[j]] = (float)(matrix[position][WeightMatrix.letters[j]] / sum);
			            }
			        }
			        // log-odds
			        for (int pos = 0; pos < matrix.length; pos++) {
			            for (int j = 0; j < WeightMatrix.letters.length; j++) {
			                matrix[pos][WeightMatrix.letters[j]] = (float)Math.log(Math.max(matrix[pos][WeightMatrix.letters[j]], .001) / 
			                		(WeightMatrix.letters[j]=='G'||WeightMatrix.letters[j]=='C'?gc/2:(1-gc)/2));
			            }
			        } 
				}
		        return wms;
			}
		}
		catch (IOException e){
			return null;
		}
    }
   
	/**
	 *  Scan the sequence for best match to the weight matrix
	 *  @return  Pair of values, the lower coordinate of highest scoring PWM hit and the score
	 *  The position will be negative if the match is on '-' strand    
	 */
	public static Pair<Integer, Double> scanPWM(String sequence, int wmLen, WeightMatrixScorer scorer){
		if (sequence==null||sequence.length()<wmLen-1){
			return new Pair<Integer, Double>(-1, Double.NEGATIVE_INFINITY);
		}
		WeightMatrixScoreProfile profiler = scorer.execute(sequence);
		double maxSeqScore = Double.NEGATIVE_INFINITY;
		int maxScoringShift = 0;
		char maxScoringStrand = '+';
		for (int i=0;i<profiler.length();i++){
			double score = profiler.getMaxScore(i);
			if (maxSeqScore<score || (maxSeqScore==score && maxScoringStrand=='-')){	// equal score, prefer on '+' strand
				maxSeqScore = score;
				maxScoringShift = i;
				maxScoringStrand = profiler.getMaxStrand(i);
			}
		}
	
		if (maxScoringStrand =='-'){
			maxScoringShift = -maxScoringShift;		
		}
		return new Pair<Integer, Double>(maxScoringShift, maxSeqScore);
	}

	/**
	 *  Scan the sequence (and reverseComplement) to find all matches to the weight matrix<br>
	 *  Note: the definition of motif position here is different from scanPWM() method<br>
	 *  position represent the match base position of the middle of the motif
	 *  
	 *  @return  List of positions (middle of motif match) that pass the threshold. <br>
	 *  The position will be negative if the match is on the reverseComplement strand     
	 */
	public static ArrayList<Integer> getAllPWMHit(String sequence, int wmLen, WeightMatrixScorer scorer, double threshold){
		ArrayList<Integer> pos = new ArrayList<Integer>();
		if (sequence==null||sequence.length()<wmLen-1){
			return pos;
		}
		WeightMatrixScoreProfile profiler = scorer.execute(sequence);
		for (int i=0;i<profiler.length();i++){
			double score = profiler.getMaxScore(i);
			if (score >= threshold){
				char maxScoreStrand = profiler.getMaxStrand_both(i);
				switch(maxScoreStrand){
				case '+':
					pos.add(i+wmLen/2);
					break;
				case '-':
					pos.add( - (i+(wmLen-1-wmLen/2)) );			// to ensure that the matching base is the same as the forward motif/strand
					break;
				case '=':	// if palindromic motif, may match on the same position in both strands
					pos.add(i+wmLen/2);
					pos.add( - (i+(wmLen-1-wmLen/2)) );	
				}
			}
		}
		return pos;
	}
	/**
	 *  Scan the sequence (and reverseComplement) to find all exact matches to the kmer<br>  
	 *  @return  List of middle positions of k-mer matches. <br>
	 *  The position represents the position of the middle of the motif match in the original sequence, <br>
	 *  mid_index = ceiling((motif length)/2), zero-based index<br>
	 *  mid_index_rc = motif length-1 - ceiling((motif length)/2), to ensure that the matching base is the same as the forward motif/strand<br>
	 *  negative position: match is on the reverseComplement strand, or, match of the RC motif on the original strand  
	 */
	public static ArrayList<Integer> getAllKmerHit(String sequence, String kmer){
		ArrayList<Integer> pos = new ArrayList<Integer>();
		if (sequence==null||kmer==null){
			return pos;
		}
		int offset_mid = kmer.length()/2;		// adjust the position from the start of the k-mer to the middle
		ArrayList<Integer> starts= new ArrayList<Integer>();
		starts = StringUtils.findAllOccurences(sequence, kmer);
		for (int p:starts)
			pos.add(p+offset_mid);
		// reverse compliment
		starts = StringUtils.findAllOccurences(sequence, SequenceUtils.reverseComplement(kmer));
		for (int p:starts)
			pos.add( - (p+kmer.length()-1-offset_mid) ); // to ensure that the matching base is the same as the forward motif/strand
		return pos;
	}	
	/**
	 *  Scan the sequence using weight matrix, outwards from the given point, until a match pass the threshold<br>
	 *  return  Pair of values, the start position of nearest PWM hit and the score<br>
	 *  The position will be negative if the match is on '-' strand    <br>
	 *  If no match pass the threshold, return -999 as position. The caller need to check for this.<br>
	 */
	public static  Pair<Integer, Double> scanPWMoutwards(String sequence, WeightMatrix wm, WeightMatrixScorer scorer, int middle, double threshold){
		if (sequence==null||sequence.length()<wm.length()-1){
			return new Pair<Integer, Double>(-999,-1.0);
		}
		WeightMatrixScoreProfile profiler = scorer.execute(sequence);
		double goodScore = 0;
		int goodScoreShift = -999;
		char goodScoreStrand = '+';
		int maxRange = Math.max(middle, profiler.length()-middle);
		for (int i=0;i<maxRange;i++){
			int idx = middle+i;
			if (idx<0 || idx>=sequence.length()-wm.length())
				continue;
			double score = profiler.getMaxScore(idx);
			if (score>=threshold){
				goodScore = score;
				goodScoreShift = idx;
				goodScoreStrand = profiler.getMaxStrand(idx);
				break;
			}
			idx = middle-i;
			if (idx<0 || idx>=sequence.length()-wm.length())
				continue;
			score = profiler.getMaxScore(idx);
			if (score>=threshold){
				goodScore = score;
				goodScoreShift = idx;
				goodScoreStrand = profiler.getMaxStrand(idx);
				break;
			}
		}
	
		if (goodScoreStrand =='-'){
			goodScoreShift = -goodScoreShift;		
		}
		return new Pair<Integer, Double>(goodScoreShift, goodScore);
	}
	
	public static void printMotifLogo(WeightMatrix wm, File f, int pixheight){
		int pixwidth = (pixheight-WeightMatrixPainter.Y_MARGIN*3-WeightMatrixPainter.YLABEL_SIZE) * wm.length() /2 +WeightMatrixPainter.X_MARGIN*2;
		System.setProperty("java.awt.headless", "true");
		BufferedImage im = new BufferedImage(pixwidth, pixheight,BufferedImage.TYPE_INT_ARGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        WeightMatrixPainter wmp = new WeightMatrixPainter();
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,pixwidth, pixheight);
        wmp.paint(wm,g2,0,0,pixwidth,pixheight);
        try {
            ImageIO.write(im,"png",f);
        }  catch (IOException ex) {
            ex.printStackTrace();
        }
	}
	
	/**
	 * Visualize sequences as color pixels
	 * @param seqs, raw sequences or FASTA sequences
	 * @param width, width of each base, in pixel
	 * @param height, height of each base, in pixel
	 * @param f, output file
	 */
	public static void visualizeSequences(String[] seqs, int width, int height, File f){
		if (seqs.length==0)
			return;
		
		int pixheight = 0;
		int maxLen = 0;
		for (String s:seqs){
        	if (s.length()!=0 && s.charAt(0)!='>')	{		// ignore header line of FASTA file
        		pixheight += height;
        		if (maxLen < s.length())
        			maxLen = s.length();
        	}
		}
		int pixwidth = maxLen*width;
		
		System.setProperty("java.awt.headless", "true");
		BufferedImage im = new BufferedImage(pixwidth, pixheight,BufferedImage.TYPE_INT_ARGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,pixwidth, pixheight);
        
        int count = 0;
        for (String s:seqs){
        	if (s.charAt(0)=='>')			// ignore header line of FASTA file
        		continue;
        	char[] letters = s.toCharArray();
        	for (int j=0;j<letters.length;j++){
        		switch(letters[j]){
        		case 'A':
        		case 'a':
        			g.setColor(Color.GREEN);
        			break;
        		case 'C':
        		case 'c':
                    g.setColor(Color.BLUE);
        			break;
        		case 'G':
        		case 'g':
                    g.setColor(Color.ORANGE);
        			break;
        		case 'T':
        		case 't':
                    g.setColor(Color.RED);
        			break;
        		case '-':
                    g.setColor(Color.WHITE);
        			break;
                default:
                	g.setColor(Color.GRAY);
        		}
                g.fillRect(j*width, count*height, width, height);
        	}
            count++;
        }
        try {
            ImageIO.write(im,"png",f);
        }  catch (IOException ex) {
            ex.printStackTrace();
        }
	}
	
	public static int mismatch(String ref, String seq){
	  if (ref.length()!=seq.length())
		  return -1;
	  int mismatch = 0;
	  for (int i=0;i<ref.length();i++){
		  if (ref.charAt(i)!=seq.charAt(i))
			  mismatch ++;
	  }
	  return mismatch;
	}
	public static void copyFile(String srFile, String dtFile){
		try{
		  File f1 = new File(srFile);
		  File f2 = new File(dtFile);
		  InputStream in = new FileInputStream(f1);
		  OutputStream out = new FileOutputStream(f2);
		
		  byte[] buf = new byte[1024];
		  int len;
		  while ((len = in.read(buf)) > 0){
			  out.write(buf, 0, len);
		  }
		  in.close();
		  out.close();
		}
		catch(FileNotFoundException ex){
		  System.err.println(ex.getMessage() + " in the specified directory.");
		  System.exit(0);
		}
		catch(IOException e){
		  System.out.println(e.getMessage());  
		}
	}
	
	/**
	 * Get a list of points that are within the window of the anchor point<br>
	 * Assuming the sites list is sorted
	 * @param sites	a list of sorted points
	 * @param anchor the anchor point
	 * @param win the window size
	 * @return
	 */
	static public ArrayList<Point> getPointsWithinWindow(ArrayList<Point> sites, Point anchor, int win){
		ArrayList<Point> results = new ArrayList<Point>();
		Region r = anchor.expand(win);
		Point start = r.startPoint();
		Point end = r.endPoint();
		int startIndex = -1;
		int endIndex = -1;
		int i = Collections.binarySearch(sites, start);
		if (i<0)
			startIndex=-i-1;		// -index-1, the insertion point
		else
			startIndex = i;
		i = Collections.binarySearch(sites, end);
		if (i<0)
			endIndex=-i-2;			// -index-1-1, the point before the insertion point
		else
			endIndex = i;
		if (startIndex<=endIndex){
			for (int j=startIndex;j<=endIndex;j++){
				results.add(sites.get(j));
			}
		}
		return results;
	}
	
	/**
	 * Count the hit count of k-mers in the sequences<br>
	 * - only count repeated kmers once in one sequence, i.e. hit count<br>
	 * - a kmer and its reverse compliment are counted as one, same for palidromic kmer 
	 * @return
	 */
	public static HashMap<String, Integer> countKmers(int k, String seqs[]){
		// expected count of kmer = total possible unique occurences of kmer in sequence / total possible kmer sequence permutation
		long tic = System.currentTimeMillis();
		
		HashMap<String, HashSet<Integer>> kmerstr2seqs = new HashMap<String, HashSet<Integer>>();
		for (int seqId=0;seqId<seqs.length;seqId++){
			String seq = seqs[seqId].toUpperCase();
			int numPos = seq.length()-k+1;
			HashSet<String> uniqueKmers = new HashSet<String>();			// only count repeated kmer once in a sequence
			 
			for (int i=0;i<numPos;i++){
				if ((i+k)>seq.length()) // endIndex of substring is exclusive
					break;
				String kstring = seq.substring(i, i+k);
				if (kstring.contains("N"))									// ignore 'N', converted from repeat when loading the sequences
					continue;
				uniqueKmers.add(kstring);
			}
			for (String s: uniqueKmers){
				if (!kmerstr2seqs.containsKey(s)){
					 kmerstr2seqs.put(s, new HashSet<Integer>());
				}
				kmerstr2seqs.get(s).add(seqId);
			}
		}
		
		// Merge kmer and its reverse compliment (RC)	
		ArrayList<String> kmerStrings = new ArrayList<String>();
		kmerStrings.addAll(kmerstr2seqs.keySet());
		
		// create kmers from its and RC's counts
		for (String key:kmerStrings){
			if (!kmerstr2seqs.containsKey(key))		// this kmer has been removed, represented by RC
				continue;
			// consolidate kmer and its reverseComplment kmer
			String key_rc = SequenceUtils.reverseComplement(key);				
			if (!key_rc.equals(key)){	// if it is not reverse compliment itself
				if (kmerstr2seqs.containsKey(key_rc)){
					int kCount = kmerstr2seqs.get(key).size();
					int rcCount = kmerstr2seqs.get(key_rc).size();
					String winner = kCount>=rcCount?key:key_rc;
					String loser = kCount>=rcCount?key_rc:key;
					kmerstr2seqs.get(winner).addAll(kmerstr2seqs.get(loser));	// winner take all
					kmerstr2seqs.remove(loser);					// remove the loser kmer because it is represented by its RC
				}
			}
		}

		HashMap<String, Integer> kmerCounts = new HashMap<String, Integer>();
		for (String key:kmerstr2seqs.keySet()){	
			kmerCounts.put(key, kmerstr2seqs.get(key).size());
		}
		kmerstr2seqs=null;
		System.gc();
		
		return kmerCounts;				
	}
	public static void printGenomeInfo (String[] args) {
		Genome genome = null;
		try {
	    	Pair<Organism, Genome> pair = Args.parseGenome(args);
	    	if(pair==null){
	    	  System.err.println("No genome provided; provide a Gifford lab DB genome name");
	    	  System.exit(1);
	    	}else{
	    		genome = pair.cdr();
	    	}
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
		Map<String, Integer> map = genome.getChromLengthMap();
		for (String chr: map.keySet())
			System.out.println(chr+"\t"+map.get(chr));
	}
	
	public static void main0(String[] args){
		System.out.println(findKey(new double[]{0,1,1,1,2,4,6}, 7));
	}
	
	public static void main1(String[] args){
		// --seqFile Y:\Tools\GPS\Multi_TFs\ESTF-Jan10\Sox2_Klf4_Tcfcp2I1_sites_MSA.fasta.txt
//	    String inName = Args.parseString(args, "seqFile", null);
//	    ArrayList<String> strs = readTextFile(inName);
//	    if (strs.isEmpty())
//	    	return;
//	    String[] ss = new String[strs.size()];
//	    strs.toArray(ss);
//	    int width = Args.parseInteger(args, "w", 5);
//	    int height = Args.parseInteger(args, "h", 3);
//	    visualizeSequences(ss, width, height, new File(inName+".png"));
	    
	}
	
    // --species "Mus musculus;mm8" --motif "CTCF" --version "090828" --windowSize 100 --motifThreshold 11.52
    public static void main3(String args[]){
		// load motif
    	Genome genome;
    	Organism org=null;
    	WeightMatrix motif = null;
    	ArgParser ap = new ArgParser(args);
		Set<String> flags = Args.parseFlags(args);		
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
	        org = pair.car();
	      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }
	    
		if (Args.parseString(args, "pfm", null)==null){
			Pair<WeightMatrix, Double> wm = CommonUtils.loadPWM(args, org.getDBID());
			motif = wm.car();
//			CommonUtils.printMotifLogo(motif, new File("test.png"), 150);	
		}
		else{
			motif = CommonUtils.loadPWM_PFM_file(Args.parseString(args, "pfm", null), Args.parseDouble(args, "gc", 0.41)); //0.41 human, 0.42 mouse
//			CommonUtils.printMotifLogo(motif, new File("test.png"), 150);
		}

		System.out.println(WeightMatrix.printMatrix(motif));
		
		System.out.println(motif.getMaxScore());
    }
    
    // testing getAllKmerHit()
    public static void main4(String args[]){
    	String s = "CATTAATTCCGTAAT";
    	String kmer = "ATTA";
    	ArrayList<Integer> pos = getAllKmerHit(s,kmer);
    	
    	System.out.println(s);
    	System.out.println(kmer);
    	for (int p:pos){
    		System.out.print(p+" ");
    	}
    	System.out.println();
    }
    
 // testing getAllPWMHit()
    public static void main(String args[]){
    	String s = "TCAGCTGTAAT";
    	List<WeightMatrix> wms = CommonUtils.loadPWMs_PFM_file("test_pwms.txt", 0.41);
    	WeightMatrix wm = wms.get(2);
    	WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
    	ArrayList<Integer> pos = getAllPWMHit(s, wm.length(), scorer, wm.getMaxScore()*0.6);
    	System.out.println(s);
    	System.out.println(WeightMatrix.printMatrixLetters(wm));
    	for (int p:pos){
    		System.out.print(p+" ");
    	}
    	System.out.println();
    }
}
