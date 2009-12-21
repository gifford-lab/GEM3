package edu.mit.csail.cgs.ewok.verbs.io;

import java.io.*;
import java.util.*;
import org.apache.log4j.Logger;

import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.CGSException;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.alignments.BLATAlignmentHitRegion;
import edu.mit.csail.cgs.datasets.alignments.PSLHitRegion;
import edu.mit.csail.cgs.datasets.alignments.MAQHitRegion;
import edu.mit.csail.cgs.datasets.alignments.parsing.BLAST8HitRegion;

public class BLATResultFileInputStream implements Iterator<BLATAlignmentHitRegion>, Closeable {

	static Logger logger = Logger.getLogger(BLATResultFileInputStream.class);

	public enum BLATResultFileFormat {
		UNKNOWN, PSL, PSLNOHEADER, BLAST8, MAQ
	}
	
	private static final int BLAST8_RESULT_NUM_TOKENS = 12;
	private static final int PSL_RESULT_NUM_TOKENS = 21;
    private static final int MAQ_RESULT_NUM_TOKENS = 16;
	
	private BufferedReader br;	
	private Genome genome;
	private BLATResultFileFormat fileFormat = null;
	private String pendingHitLine = null;
    private String source;

	/**
	 * 
	 * @param filename
	 * @param organism
	 * @param genome
	 * @throws CGSException
	 */
	public BLATResultFileInputStream(String filename, Genome genome) throws CGSException {
        try {
            init(new BufferedReader(new InputStreamReader(new FileInputStream(filename))), genome);
        } catch (FileNotFoundException e) {
            throw new CGSException(e);
        }
        source = filename;
	}

	/**
	 * 
	 * @param f
	 * @param organism
	 * @param genome
	 * @throws CGSException
	 */
	public BLATResultFileInputStream(File f, Genome genome) throws CGSException {		
        try {
            init(new BufferedReader(new InputStreamReader(new FileInputStream(f))), genome);
        } catch (FileNotFoundException e) {
            throw new CGSException(e);
        }
        source = f.toString();
    }

	public BLATResultFileInputStream(BufferedReader r, Genome genome) throws CGSException {		
        init(r,genome);
        source = r.toString();
    }
        
    private void init(BufferedReader reader, Genome genome) throws CGSException {
		this.genome = genome;
		try {
			br = reader;
			pendingHitLine = null;
			String line;
			line = br.readLine();
			if (line != null) {
				this.determineFileFormat(line);
				switch (fileFormat) {
				case PSL: this.readThroughPSLHeader();
                    line = br.readLine();
                    pendingHitLine = line;
                    break;
                case PSLNOHEADER:
                    pendingHitLine = line;
                    break;
				case BLAST8: pendingHitLine = line;
                    break;
                case MAQ:
                    pendingHitLine=line;
                    break;
				case UNKNOWN: pendingHitLine = null;
                    break;
				}
			}
		}
		catch (IOException ioex) {
			throw new CGSException(ioex);
		}
	}


	/**
	 * 
	 */
	public boolean hasNext() {
		return (pendingHitLine != null);
	}
    public void remove() throws UnsupportedOperationException {
        throw new UnsupportedOperationException("Can't remove from a BLATResultFileInputStream");
    }

	/**
	 * 
	 */
	public BLATAlignmentHitRegion next() {
		BLATAlignmentHitRegion hit = null;
		try {
			switch (fileFormat) {
			case PSL: 
            case PSLNOHEADER:
                hit = parsePSLResult(pendingHitLine);
                break;
			case BLAST8: 
                hit = parseBLAST8Result(pendingHitLine);
                break;
            case MAQ:
                hit = parseMAQResult(pendingHitLine);
                break;
			}
		}
		catch (CGSException cgsex) {
			//TODO modify the method in the interface to allow an exception to be thrown
			logger.error(cgsex,cgsex);
			hit = null;
		}
        /* read the next line, even if hit = null.  This
           lets you skip misformarted parts of the data and
           keep going
        */
        try {
            pendingHitLine = br.readLine();
            if (pendingHitLine == null) {
                br.close();
            }
        }
        catch (IOException ioex) {
            logger.error(ioex);
        }
		return hit;
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.mit.csail.cgs.hyperdrive.utils.Stream#getDescriptor()
	 */
	public String getDescriptor() {
		return "BLAT Results(" + source + ")";
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.mit.csail.cgs.utils.Closeable#close()
	 */
	public void close() {
		if (isClosed()) {
			return;
		}
		try {
			br.close();
		} catch (IOException ioex) {
			logger.error("Error closing buffered reader", ioex);

			// TODO throw a generic (e.g. CGS) exception

		}
		br = null;
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
	 */
	public boolean isClosed() {
		return br == null;
	}


	/**
	 * 
	 * @param firstLine
	 */
	private void determineFileFormat(String firstLine) {
        int fields = firstLine.split("\t").length;
		if (firstLine.startsWith("psLayout version 3")) {
			fileFormat = BLATResultFileFormat.PSL;             
		} else if (fields == PSL_RESULT_NUM_TOKENS) {
            fileFormat = BLATResultFileFormat.PSLNOHEADER;
        } else {
			try {
				BLATAlignmentHitRegion testHit = parseBLAST8Result(firstLine);
				fileFormat = BLATResultFileFormat.BLAST8;				
			}
			catch (CGSException cgsex) {
                try {
                    BLATAlignmentHitRegion testHit = parseMAQResult(firstLine);
                    fileFormat = BLATResultFileFormat.MAQ;				
                } catch (CGSException cgsex2) {
                    fileFormat = BLATResultFileFormat.UNKNOWN;
                }
			}
		}
	}
	
	
	/**
	 * Attempts to read four lines - the remainder of the PSL header after the
	 * first header line
	 * @throws IOException if an error occurs while reading a line
	 */
	private void readThroughPSLHeader() throws IOException {
		String line = br.readLine();
		int numRead = 1;
		while ((numRead < 4) && (line != null)) {
			br.readLine();
			numRead++;
		}
	}
	
	
	/**
	 * 
	 * @param resultLine
	 * @return
	 * @throws CGSException
	 */
	private BLATAlignmentHitRegion parseBLAST8Result(String resultLine) throws CGSException {
		return this.parseBLAST8Result(resultLine, null);
	}
	
	
	/**
	 * 
	 * @param resultLine
	 * @param queryIDPrefix
	 * @return
	 * @throws CGSException
	 */
	private BLATAlignmentHitRegion parseBLAST8Result(String resultLine, String queryIDPrefix) throws CGSException {
		//TODO make the queryIDPrefix a regex pattern
		String[] tokens = resultLine.split("\\s+");
		if (tokens.length == BLAST8_RESULT_NUM_TOKENS) {
			try { 
				String queryID = tokens[0];
				if ((queryIDPrefix != null) && !queryID.startsWith(queryIDPrefix)) {
					throw new CGSException("Query ID does not start with specified query ID prefix.");
				}
				String chrom = tokens[1];
				int subjectStart = Integer.valueOf(tokens[8]);
				int subjectEnd = Integer.valueOf(tokens[9]);
				char strand;
				if (subjectStart < subjectEnd) {
					strand = '+';
				}
				else {
					strand = '-';
                    subjectStart = Integer.valueOf(tokens[9]);
                    subjectEnd = Integer.valueOf(tokens[8]);
				}
				
				double percentIdentity = Double.valueOf(tokens[2]);
				int alignmentLength = Integer.valueOf(tokens[3]);
				int numMismatches = Integer.valueOf(tokens[4]);
				int numGapOpenings = Integer.valueOf(tokens[5]);
				int queryStart = Integer.valueOf(tokens[6]) - 1;
				int queryEnd = Integer.valueOf(tokens[7]) - 1;
				
				double eValue = Double.valueOf(tokens[10]);
				double bitScore = Double.valueOf(tokens[11]);
	
				
				BLATAlignmentHitRegion hit = 
					new BLAST8HitRegion(genome, chrom, subjectStart, 
                                        subjectEnd, queryID, strand, 
                                        percentIdentity, alignmentLength, numMismatches, 
                                        numGapOpenings, queryStart, queryEnd,  eValue, bitScore);

			
				return hit;
			}
			catch (NumberFormatException nfex) {
				throw new CGSException(nfex);
			}
		}
		else {
            System.err.println("BAD LINE " + resultLine);
			throw new CGSException("Result line has " + tokens.length + " tokens instead of " + BLAST8_RESULT_NUM_TOKENS + ".");
		}
	}
	

    private BLATAlignmentHitRegion parseMAQResult(String resultLine) throws CGSException {
        String[] tokens = resultLine.split("\\t");
        if (tokens.length == MAQ_RESULT_NUM_TOKENS) {
            try {
                return new MAQHitRegion(genome,
                                        tokens[0],
                                        tokens[1],
                                        Integer.parseInt(tokens[2]),
                                        tokens[3].charAt(0),
                                        Integer.parseInt(tokens[4]),
                                        Integer.parseInt(tokens[5]) == 1,
                                        Integer.parseInt(tokens[6]),
                                        Integer.parseInt(tokens[7]),
                                        Integer.parseInt(tokens[8]),
                                        Integer.parseInt(tokens[9]),
                                        Integer.parseInt(tokens[13]));
            } catch (NumberFormatException nfex) {
				throw new CGSException(nfex);
			}
        } else {
            throw new CGSException("Result line has " + tokens.length + " tokens instead of " + MAQ_RESULT_NUM_TOKENS + ": " + resultLine);
        }
    }
	
	/**
	 * 
	 * @param resultLine
	 * @return
	 * @throws CGSException
	 */
	private BLATAlignmentHitRegion parsePSLResult(String resultLine) throws CGSException {
		return this.parsePSLResult(resultLine, null);
	}
	
	
	/**
	 * 
	 * @param resultLine
	 * @return
	 */
	private BLATAlignmentHitRegion parsePSLResult(String resultLine, String queryNamePrefix) throws CGSException {
		//TODO make the queryNamePrefix a regex pattern
		String[] tokens = resultLine.split("\\t");
		if (tokens.length == PSL_RESULT_NUM_TOKENS) {
			try { 
				int matches = Integer.valueOf(tokens[0]); //number of bases that match that aren't repeats
				int numMismatches = Integer.valueOf(tokens[1]); //number of bases that don't match
				int repMatches = Integer.valueOf(tokens[2]); //number of bases that match but are part of repeats
				int nCount = Integer.valueOf(tokens[3]); //number of 'N' bases
				int qNumInsert = Integer.valueOf(tokens[4]); //Number of inserts in query
				int qBaseInsert = Integer.valueOf(tokens[5]); //Number of bases inserted in query
				int tNumInsert = Integer.valueOf(tokens[6]); //Number of inserts in target
				int tBaseInsert = Integer.valueOf(tokens[7]); //Number of bases inserted in target
				String strandString = tokens[8]; //+ or - for query strand, optionally followed by + or - for target strand
				String qName = tokens[9]; //Query sequence name
				int qSize = Integer.valueOf(tokens[10]); //Query Sequence size
			    int qStart = Integer.valueOf(tokens[11]); //Alignment start position in query
			    int qEnd = Integer.valueOf(tokens[12]) - 1; //Alignment end position in query
			    String chrom = tokens[13]; //Target sequence name
			    int tSize = Integer.valueOf(tokens[14]); //Target sequence size
			    int tStart = Integer.valueOf(tokens[15]); //Alignment start position in target
			    int tEnd = Integer.valueOf(tokens[16]); //Alignment end position in target
			    int blockCount = Integer.valueOf(tokens[17]); //Number of blocks in alignment. A block contains no gaps.
			    String blockSizesString = tokens[18]; //Size of each block in a comma separated list
			    String qStartsString = tokens[19]; //Start of each block in query in a comma separated list
			    String tStartsString = tokens[20]; //Start of each block in target in a comma separated list
				
			    char strand = strandString.charAt(0);
				if (strandString.length() == 2) {
					strand = strandString.charAt(1);
				}
				
				if ((queryNamePrefix != null) && !qName.startsWith(queryNamePrefix)) {
					throw new CGSException("Query Name does not start with specified query ID prefix.");
				}
				
				//parse the block sizes into a list
				String[] blockSizeTokens = blockSizesString.split(",");
				List<Integer> blockSizes = new Vector<Integer>();
				for (int i = 0; i < blockSizeTokens.length; i++) {
					blockSizes.add(Integer.valueOf(blockSizeTokens[i]));
				}
				
				//parse the query start sites into a list
				String[] qStartTokens = qStartsString.split(",");
				List<Integer> qStarts = new Vector<Integer>();
				for (int i = 0; i < qStartTokens.length; i++) {
					qStarts.add(Integer.valueOf(qStartTokens[i]));
				}

				//parse the target start sites into a list
				String[] tStartTokens = tStartsString.split(",");
				List<Integer> tStarts = new Vector<Integer>();
				for (int i = 0; i < tStartTokens.length; i++) {
					tStarts.add(Integer.valueOf(tStartTokens[i]));
				}
				
				//compute some values that aren't directly specified in PSL format
				double percentIdentity = (double)(matches + repMatches) / (double)(matches + repMatches + numMismatches);
				int alignmentLength = tEnd - tStart;
				int numGapOpenings = qNumInsert + tNumInsert;
				
				
				//create the hit object
				BLATAlignmentHitRegion hit = 
					new PSLHitRegion(genome, chrom, tStart, tEnd, qName, strand, 
							percentIdentity, alignmentLength, numMismatches, 
							numGapOpenings, qStart, qEnd, 
							matches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert,
							qSize, tSize, blockCount, blockSizes, qStarts, tStarts);
				
				return hit;
			}
			catch (NumberFormatException nfex) {
				throw new CGSException(nfex);
			}
		}
		else {
			throw new CGSException("Result line has " + tokens.length + " tokens instead of " + BLAST8_RESULT_NUM_TOKENS + ": " + resultLine);
		}
	}
}
