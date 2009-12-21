package edu.mit.csail.cgs.datasets.alignments.parsing;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.parsing.FASTAStream;

public class BLAST8Parser {

	static Logger logger = Logger.getLogger(BLAST8Parser.class);

	/**
	 * Parses data in the blast8 format, e.g.
	 * WICMT-SOLEXA_10:2:300:41:523    chr1    100.00  19      0       0       7       
	 * 25      110891815       110891833       3.8e-02 36.0
	 * WICMT-SOLEXA_10:2:300:41:523    chr1    100.00  6       0       0       1
	 * 6       110891604       110891609       1.1e+06 11.0
	 *
	 * @param filename name of the file containing the data
	 * @return a List of hit objects
	 */
	public static List<BLAST8Hit> parseBLAST8Output(String filename) {
		Vector<BLAST8Hit> results = new Vector<BLAST8Hit>();

		FileReader in = null;
		BufferedReader bin = null;
		try {
			in = new FileReader(filename);
			bin = new BufferedReader(in);

			String hitString = bin.readLine();
			while (hitString != null) {
				BLAST8Hit hit = BLAST8Parser.parseLine(hitString, 0);
				results.add(hit);
				hitString = bin.readLine();
			}
		}
		catch(IOException ioex) {
			logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				logger.error("Error closing buffered reader", ioex2);
			}			
		}
		return results;
	}


	/**
	 * 
	 * @param filename
	 */
	public static void outputAlignmentLengthStats(String filename) {
		TreeMap<Integer, Long> lengthHist = new TreeMap<Integer, Long>();

		FileReader in = null;
		BufferedReader bin = null;
		try {
			in = new FileReader(filename);
			bin = new BufferedReader(in);

			String hitString = bin.readLine();
			int i = 1;
			while (hitString != null) {
				BLAST8Hit hit = BLAST8Parser.parseLine(hitString, i);
				if (hit != null) {
					if (!lengthHist.containsKey(hit.alignLength)) {
						lengthHist.put(hit.alignLength, (long)0);
					}
					lengthHist.put(hit.alignLength, (lengthHist.get(hit.alignLength) + (long)1));
				}

				hitString = bin.readLine();				
				i++;
			}
			for (Map.Entry<Integer, Long> entry : lengthHist.entrySet()) {
				System.out.println(entry.getKey() + ", " + entry.getValue());
			}		
		}
		catch(IOException ioex) {
			logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				logger.error("Error closing buffered reader", ioex2);
			}			
		}
	}

	/**
	 * Parse a single line of text into a hit object
	 * @param blast8Line a line of text representing a hit
	 * @return a hit object containing the data from the specified line
	 */
	private static BLAST8Hit parseLine(String blast8Line, int lineNumber) {
		BLAST8Hit hit = new BLAST8Hit();
		String[] tokens = blast8Line.split("\\s");
		if (tokens.length == 12) {
			try { 
				hit.queryID = tokens[0];
				if (!hit.queryID.startsWith("WICMT")) {
					logger.error("Parse error on line " + lineNumber + " query has invalid query ID " + hit.queryID + ".");
					return null;
				}
				hit.subjectID = tokens[1];
				if (!hit.subjectID.startsWith("chr")) {
					logger.error("Parse error on line " + lineNumber + " query has invalid subject ID " + hit.queryID + ".");
					return null;
				}
				hit.percentIdentity = Double.valueOf(tokens[2]);
				hit.alignLength = Integer.valueOf(tokens[3]);
				hit.mismatches = Integer.valueOf(tokens[4]);
				hit.gapOpenings = Integer.valueOf(tokens[5]);
				hit.queryStart = Integer.valueOf(tokens[6]);
				hit.queryEnd = Integer.valueOf(tokens[7]);
				hit.subjectStart = Integer.valueOf(tokens[8]);
				hit.subjectEnd = Integer.valueOf(tokens[9]);
				hit.eValue = Double.valueOf(tokens[10]);
				hit.bitScore = Double.valueOf(tokens[11]);
			}
			catch (Exception ex) {
				logger.error("Parse error on line " + lineNumber + ".", ex);
				return null;
			}
		}
		else {
			logger.error("Line " + lineNumber + " has " + tokens.length + " tokens.");
			return null;
		}
		return hit;
	}


	public static void outputCorruptedResultSequences(String[] fastaFilenames, String[] blatFilenames) {
		for (int i = 0; i < blatFilenames.length; i++) {
			Vector<Pair<String, String>> corruptSequences = new Vector<Pair<String, String>>();

			logger.error("Parsing FASTA File: " + fastaFilenames[i]);
			logger.error("Parsing blat file: " + blatFilenames[i]);
			File fastaFile = new File(fastaFilenames[i]);
			FASTAStream fastaStream = null;
			FileReader in = null;
			BufferedReader bin = null;
			try {
				//open fasta file
				fastaStream = new FASTAStream(fastaFile);
				Pair<String, String> currentFASTASequence = null;
				if (fastaStream.hasNext()) {
					currentFASTASequence = fastaStream.next();
				}
				else {
					System.out.println("Can't get first sequence ID from FASTA file");
					System.exit(0);
				}

				//open blat file
				in = new FileReader(blatFilenames[i]);
				bin = new BufferedReader(in);

				String hitString = bin.readLine();
				int lineNumber = 1;
				while (hitString != null) {
					BLAST8Hit hit = BLAST8Parser.parseLine(hitString, lineNumber);
					if (hit != null) {
						if (!hit.queryID.equals(currentFASTASequence.getFirst())) {
						    boolean done = false;
							while (fastaStream.hasNext() && !done) {
								currentFASTASequence = fastaStream.next();
								if (hit.queryID.equals(currentFASTASequence.getFirst())) {
								    done = true;
								}
							}
							if (!hit.queryID.equals(currentFASTASequence.getFirst())) {
									System.out.println("ERROR: FASTA Sequence ID does not match BLAT Query ID!\nFASTA Sequence ID: " 
											+ currentFASTASequence.getFirst() + "\nBLAT Query ID: " + hit.queryID);
									System.exit(0);
							}
						}
					}
					else {
						//keep reading lines until a good line is reached
					    logger.error("Corrupt sequence after hit with ID: " + currentFASTASequence.getFirst());
					    if (!corruptSequences.contains(currentFASTASequence)) {
						corruptSequences.add(currentFASTASequence);
					    }
						BLAST8Hit goodHit = null;
						while ((hitString != null) && (goodHit == null)) {
							hitString = bin.readLine();
							lineNumber++;
							goodHit = BLAST8Parser.parseLine(hitString, lineNumber);								
						}
						if (goodHit != null) {
						    logger.error("Found good hit on line: " + lineNumber + " with query ID: " + goodHit.queryID);
						    logger.error("Current FASTA Sequence ID: " + currentFASTASequence.getFirst());
							boolean done = false;
							if (goodHit.queryID.equals(currentFASTASequence.getFirst())) {
							    logger.error("Good hit has same query ID as previous good hit");
							    done = true;
							}
							while (fastaStream.hasNext() && !done) {
								Pair<String, String> nextSequence = fastaStream.next(); 
								corruptSequences.add(nextSequence);
								if (goodHit.queryID.equals(nextSequence.getFirst())) {
									done = true;
									logger.error("Found query ID match");
								}
								currentFASTASequence = nextSequence;
							}
						}
						else {
							//read to the end of the BLAT file
						    logger.error("Read to end of BLAT file without finding a good hit");
						    while (fastaStream.hasNext()) {
								Pair<String, String> nextSequence = fastaStream.next(); 
								corruptSequences.add(nextSequence);
							}
							break; //break out of the loop reading through the blat file
						}
					}

					hitString = bin.readLine();				
					lineNumber++;
				}

				if (corruptSequences.size() > 0) {
					PrintStream ps = null;
					FileOutputStream os = null;
					int lineLength;

					try {
						String corruptSequenceFilename = fastaFilenames[i] + "_corrupt";
						os = new FileOutputStream(corruptSequenceFilename);
						ps = new PrintStream(os);
						lineLength = 80;

						for (Pair<String, String> query : corruptSequences) {
							ps.println(">" + query.getFirst());
							String s = query.getLast();
							int start = 0;
							int len = s.length();
							while (start < len) {
								int end = start + lineLength;
								if (end > len) {
									end = len;
								}
								ps.println(s.substring(start,end));
								start+= lineLength;
							}							 
						}
					}
					catch (IOException ioex) {
						logger.error("Error writing new file", ioex);
					}
					finally {
						try {
							ps.close();
							if (os != null) {
								os.close();                
							}
							ps.close();
						}
						catch(IOException ioex2) {
							//nothing left to do here, just log the error
							logger.error("Error closing buffered reader", ioex2);
						}		
					}
				}
			}
			catch(IOException ioex) {
				logger.error("Error parsing file", ioex);
			}
			finally {
				try {
					if (bin != null) {
						bin.close();
					}
				}
				catch(IOException ioex2) {
					//nothing left to do here, just log the error
					logger.error("Error closing buffered reader", ioex2);
				}		
				if (fastaStream != null) {
					fastaStream.close();
				}
			}
		}
	}


	public static void main(String[] args) {	
		String[] fastaFilenames = {"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_aa",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ab",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ac",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ad",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ae",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_af", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ag",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ah",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ai",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_aj",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ak",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_al",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_am",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_an",
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ao", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ap", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_aq", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ar", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_as", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_at", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_au", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_av", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_aw", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ax", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ay", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_az", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ba", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bb", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bc", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bd", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_be", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bf", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bg", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bh", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bi", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bj", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bk", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bl", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bm", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bn", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bo", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bp", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bq", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_br", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bs", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bt", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bu", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bv", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bw", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bx", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_by", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_bz", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ca", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_cb", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_cc", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_cd", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ce", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_cf", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_cg", 
				"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ch", 
		"/scratch/mes_oct4_solexa/blat/split60/mes_oct4_solexa.fasta_ci"};

		
		//BLAST8Parser.parseBLAST8Output("G:\\projects\\chip_seq\\tmp\\mes_oct4_solexa.fasta_ci.against_mm8.blat");
		String[] blatFilenames = {"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_aa.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ab.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ac.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ad.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ae.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_af.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ag.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ah.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ai.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_aj.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ak.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_al.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_am.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_an.against_mm8.blat",
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ao.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ap.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_aq.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ar.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_as.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_at.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_au.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_av.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_aw.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ax.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ay.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_az.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ba.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bb.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bc.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bd.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_be.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bf.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bg.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bh.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bi.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bj.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bk.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bl.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bm.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bn.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bo.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bp.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bq.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_br.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bs.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bt.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bu.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bv.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bw.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bx.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_by.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_bz.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ca.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_cb.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_cc.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_cd.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ce.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_cf.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_cg.against_mm8.blat", 
				"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ch.against_mm8.blat", 
		"/scratch/mes_oct4_solexa/blat/mes_oct4_solexa.fasta_ci.against_mm8.blat"};

		BLAST8Parser.outputCorruptedResultSequences(fastaFilenames, blatFilenames);
		
//		for (int i = 0; i < blatFilenames.length; i++) {
//			System.out.println("Parsing " + blatFilenames[i]);
//			BLAST8Parser.outputAlignmentLengthStats(blatFilenames[i]);
//		}		
	}
}
