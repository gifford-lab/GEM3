package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class MergeKmerData {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ArrayList<String> texts = CommonUtils.readTextFile(Args.parseString(args, "design", "design.txt"));
		// count number of data columns
		int count = 0;
		for (String line:texts){
			if (line.length()==0)
				continue;
			String[] f = line.split("\t");
			String[] cols = f[2].split(",");
			count += cols.length;
		}
		
		HashMap<String, String[]> kmer2data = new HashMap<String, String[]>();
		
		StringBuilder out = new StringBuilder("Kmer\t");
		
		int col_id = 0;
		for (String line:texts){
			if (line.length()==0)
				continue;
			String[] f = line.split("\t");
			String fileName = f[0];
			String format = f[1];
			String[] cols = f[2].split(",");
			int[] ids = new int[cols.length];			// column ids to be output
			for (int i=0;i<cols.length;i++){
				ids[i]=Integer.parseInt(cols[i]);
			}
			
			int headerRows=0;
			if (format.endsWith("GEM")){
				headerRows=3;
			} 
			else if (format.endsWith("PBM")){
				headerRows=1;
			}
			
			ArrayList<String> data = CommonUtils.readTextFile(fileName);
			String header = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
			if (headerRows!=0)
				header = data.get(headerRows-1);
			String[] hf = header.split("\t");		// header field
			for (int id: ids)
				out.append(hf[id]+"\t");
			for (int i=headerRows;i<data.size();i++){
				String[] df = data.get(i).split("\t");	// data field
				String kmer = df[0];
				
				if (kmer.contains("."))			// ignore PBM gapped k-mers
					continue;
				
				if (kmer.contains("/")){
					String[] kf=kmer.split("/");
					if (kmer2data.containsKey(kf[0]))
						kmer = kf[0];
					else if (kmer2data.containsKey(kf[1]))
						kmer = kf[1];
					else	// if new k-mer
						kmer = kf[0];
				}
				else {
					String kmer_rc = SequenceUtils.reverseComplement(kmer);
					if (kmer2data.containsKey(kmer_rc))
						kmer = kmer_rc;
				}
				if (!kmer2data.containsKey(kmer))
					kmer2data.put(kmer, new String[count]);
				
				for (int j=0;j<ids.length;j++){
					String[] dataRow = kmer2data.get(kmer);
					dataRow[col_id+j]=df[ids[j]];
				}
			}
			col_id += ids.length;
		}
		out.append("\n");
		
		// output
		for (String kmer : kmer2data.keySet()){
			String[] data = kmer2data.get(kmer);
			out.append(kmer+"\t");
			for (int i=0;i<data.length;i++){
				if (data[i]==null)
					data[i]="0";
				out.append(data[i]+"\t");
			}
			out.append("\n");
		}
		CommonUtils.writeFile(Args.parseString(args, "out", "out.txt"), out.toString());

	}

}
