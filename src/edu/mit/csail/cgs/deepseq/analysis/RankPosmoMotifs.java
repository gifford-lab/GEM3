package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class RankPosmoMotifs {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ArrayList<String> text = CommonUtils.readTextFile(args[0]);
		ArrayList<Integer> hitCounts = new ArrayList<Integer>();
		for (String s: text){
			int idx = s.indexOf("Total: ");
			if (idx!=-1){
				String subs = s.substring(idx+7, s.length());
				int idxEnd = subs.indexOf(" ");
				if (idxEnd!=-1){
					hitCounts.add(Integer.parseInt(subs.substring(0, idxEnd)));
				}
			}
		}
		int[] hit_sorted = new int[hitCounts.size()];
		for (int i=0;i<hitCounts.size();i++){
			hit_sorted[i]=hitCounts.get(i);
		}
		StatUtil.findSort(hit_sorted);
		String[] posmo_motifs = new String[hit_sorted.length];
		StringBuilder motif = new StringBuilder();
		int idx_motif = -1;
		for (String s: text){
			int idx = s.indexOf("Total: ");
			if (idx!=-1){	// if header line
				if (idx_motif!=-1){	// not first line, store motif of previous block
					posmo_motifs[idx_motif]=motif.toString();
					hit_sorted[idx_motif]=-1;		// mark as used
				}
				String subs = s.substring(idx+7, s.length());
				int idxEnd = subs.indexOf(" ");
				if (idxEnd!=-1){
					 int hit = Integer.parseInt(subs.substring(0, idxEnd));
					 for (int i=0;i<hit_sorted.length;i++){
						 if (hit_sorted[i]==hit){
							 idx_motif = i;
							 motif = new StringBuilder(s).append("\n");				// start with header line
						 }							 
					 }
				}
			}
			else{		// if data line, append
				if (s.trim().equals(""))
					continue;
				motif.append(s).append("\n");
			}
		}
		if (idx_motif != -1){	// store last block
			posmo_motifs[idx_motif]=motif.toString();
			hit_sorted[idx_motif]=-1;		// mark as used
		}
		
		// output
		motif = new StringBuilder();
		for (int i = posmo_motifs.length-1;i>=0;i--){		// reverse order, print descending
			motif.append(posmo_motifs[i]);
		}
		CommonUtils.writeFile(args[0].replace(".txt", "_sorted.txt"), motif.toString());
	}

}
