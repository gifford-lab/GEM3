package edu.mit.csail.cgs.tools.sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

/**
 * This class takes input in <tt>Genbank</tt> format and converts them into output in <tt>Fasta</tt> format<br>
 * The sequences are put in windows of 60 residues by default.
 * @author gio_fou
 *
 */
public class GenbankToFasta {
	
	private static int window = 60;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		try {
			String ipfs = "/afs/csail.mit.edu/group/psrg/projects/sigma/LitData/Miura/Data/DDBJ/ESTs_s288c_sk1_31847.genbank";
			String opfs = "/Users/gio_fou/Desktop/io_tests/other_ESTs_s288c_sk1_31847.fasta";
			File ipf = new File(ipfs);
			File opf = new File(opfs);
			
			FileInputStream is = new FileInputStream(ipf);
			FileOutputStream os = new FileOutputStream(opf); 
			
			// GenbankToFasta.convertGenbankToFasta(is, os);
			
			GenbankToFasta.convertGenbankToFasta(is, os, "GENBANK", "DNA");
			
			System.out.println("THIS IS THE END!");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Converts inputs in <tt>Genbank</tt> to outputs in <tt>Fasta</tt> format<br>
	 * Except for the sequence id, it also stores other info as strain, locus, length, type etc.
	 * @param is The input file in <tt>Genbank</tt> format. E.g. <tt>System.in</tt>, <tt>new FileInputStream(File f)</tt>
	 * @param os The output file in <tt>Fasta</tt> format. E.g. <tt>System.out</tt>, <tt>new FileOutputStream(File f)</tt> 
	 */
	public static void convertGenbankToFasta(InputStream is, OutputStream os)
	{
		BufferedReader br = null;
		BufferedWriter bw = null;
		try 
		{
			InputStreamReader isr = new InputStreamReader(is);
			br = new BufferedReader(isr);
			
			OutputStreamWriter osw = new OutputStreamWriter(os);
            bw = new BufferedWriter(osw);
            
			SequenceIterator sequenceItr = SeqIOTools.readGenbank(br);
			
			int count = 0;
			while(sequenceItr.hasNext())
			{
				StringBuilder sbName = new StringBuilder();
				Sequence seq = sequenceItr.nextSequence();
				
				// get the name of the sequence
				String seqName = seq.getName();
				// get the char-sequence of the sequence
				String seqString = seq.seqString();
				
				//present sequence in steps of window (default: 60)
				String seqStringFormatted = formatSequence(seqString, window);
				
				sbName.append(seqName);
				
				String str;
				// append strain
				Iterator<Feature> featItr = seq.features();
				while( featItr.hasNext() )
				{
					Feature feat = featItr.next();
					Annotation featAnnot = feat.getAnnotation();
					if(featAnnot.containsProperty("strain"))
						if( (str = (String)featAnnot.getProperty("strain")) != null )
							sbName.append(String.format("|%s", str));
				}
				
				String[] keys = {"LOCUS", "SIZE", "TYPE", "CIRCULAR", "DIVISION", "MDAT", "SOURCE"};
				Annotation seqAnnot = seq.getAnnotation();
				for(String key : keys)
					if(seqAnnot.containsProperty(key))
						if( (str = (String)seqAnnot.getProperty(key)) != null )
							sbName.append(String.format("|%s", str));

				// Print record
				String seqRecord = String.format(">%s%n%s%n", sbName, seqStringFormatted);
				bw.write(seqRecord);
				
				// print statistics
				count++;
				if( count%1000 == 0)
					System.out.printf("Record %d has been read...%n", count);
			}					
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		catch (NoSuchElementException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
		catch (BioException e) {
			e.printStackTrace();
			System.exit(-1);
		}		
		finally
		{
			try {
				br.close();
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}// end of convertGenbankToFasta method
	
	/**
	 * Converts inputs in <tt>Genbank</tt> to outputs in <tt>Fasta</tt> format<br>
	 * Only the name of the sequence is stored on the header section.  
	 * @param is The input file in <tt>Genbank</tt> format. E.g. <tt>System.in</tt>, <tt>new FileInputStream(File f)</tt>
	 * @param os The output file in <tt>Fasta</tt> format. E.g. <tt>System.out</tt>, <tt>new FileOutputStream(File f)</tt>
	 * @param format The format of the file. Allowed formats are (case insensitive): FASTA, EMBL, GENBANK, SWISSPROT (or swiss), GENPEPT
	 * @param alpha The sustance of the file. Allowed formats are (case insensitive): DNA, AA, Protein, RNA
	 */
	public static void convertGenbankToFasta(InputStream is, OutputStream os, String format, String alpha)
	{
		if( !(format.equals("FASTA") | format.equals("EMBL") | format.equals("GENBANK") | format.equals("SWISSPROT") | format.equals("swiss") | format.equals("GENPEPT") ) )
			throw new IllegalArgumentException("Illegal value for argument format");
		
		if( !(alpha.equals("DNA") | alpha.equals("AA") | alpha.equals("Protein") | alpha.equals("RNA") ) )
			throw new IllegalArgumentException("Illegal value for argument alpha");
		
		BufferedReader br = null;
		
		try {
			InputStreamReader isr = new InputStreamReader(is);
			br = new BufferedReader(isr);
			
			SequenceIterator iter = (SequenceIterator)SeqIOTools.fileToBiojava(format, alpha, br);
		    SeqIOTools.writeFasta(os, iter);	
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		catch (NoSuchElementException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
		catch (BioException e) {
			e.printStackTrace();
			System.exit(-1);
		}		
		finally {
			try {
				br.close();
				os.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}// end of convertGenbankToFasta method

	
	/**
	 * Formats the input string in windows of size <tt>window</tt>
	 * @param seq string to be formatted
	 * @param window size of the window
	 * @return the formatted string
	 */
	public static String formatSequence(String seq, int window)
	{
		StringBuilder seqFormatted = new StringBuilder();
		int currStart = 0;
		while( currStart < seq.length())
		{
			int currEnd = Math.min(currStart+window, seq.length());
			seqFormatted.append(seq.substring(currStart, currEnd) + "\n");
			currStart = currEnd;
		}
		
		return seqFormatted.toString();
	}//end of formatSequence method
	
	/**
     * Formats the input string with the default window size of 60 residues
	 * @param seq string to be formatted
	 * @return
	 */
	public static String formatSequence(String seq)
	{
		return formatSequence(seq, window);
	}//end of formatSequence method
	
}//end of GenbankToFasta class
