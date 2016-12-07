package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Pattern;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC0;

/**
 * This motif formatter is to re-format HMS output motif to STAMP format
 * @author yuchun
 *
 */
public class MotifFormatter {
	
	public static void main(String[] args) {
		Set<String> flags = Args.parseFlags(args);
		File file = new File(args[0]);
		float[][] freqs = new float[4][];
		try{
		    BufferedReader br = new BufferedReader(new FileReader(file));
		    String line;
		    WeightMatrix matrix = null;
		    int motifCount = 0;
		    Vector<float[]> arrays = new Vector<float[]>();

		    // Read in Transfac format first
		    Organism currentSpecies = null;
		    String name = null, id = null, accession = null;
		    int lineCount = 0;
		    while ((line = br.readLine()) != null) {
		      line = line.trim();
		      if (line.length() > 0) {
		        String[] pieces = line.split("\\s+");
		        float[] fs = new float[pieces.length];
		        for(int i=0;i<fs.length;i++){
		        	fs[i]=Float.parseFloat(pieces[i]);
		        }
		        freqs[lineCount]=fs;
		      }
		      lineCount++;
		    }
		}
		catch (IOException e){
			System.out.println(file.getName()+" motif PFM file reading error!!!");
		}
		float[][] pfm = new float[freqs[0].length][freqs.length];
		for (int i=0;i<pfm.length;i++)
			for (int j=0;j<freqs.length;j++)
				pfm[i][j] = freqs[j][i];
		System.out.println(makeTRANSFAC(pfm, "DE Motif1\n"));
	}
	
	private static String makeTRANSFAC (float[][] pfm, String header){
		// make string in TRANSFAC format
		StringBuilder sb = new StringBuilder();
		sb.append(header);
		for (int p=0;p<pfm.length;p++){
			sb.append(p).append(" ");
			int maxBase = 0;
			float maxCount=0;
			for (int b=0;b<KMAC0.LETTERS.length;b++){
				sb.append(String.format("%.5f ", pfm[p][b]));
				if (maxCount<pfm[p][b]){
					maxCount=pfm[p][b];
					maxBase = b;
				}
			}
			sb.append(KMAC0.LETTERS[maxBase]).append("\n");
		}
		sb.append("XX\n\n");
		return sb.toString();
	}

}
