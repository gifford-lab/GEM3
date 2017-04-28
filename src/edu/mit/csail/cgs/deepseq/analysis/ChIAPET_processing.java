package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class ChIAPET_processing {
	// --fq1 /Users/yguo/Desktop/PPG_Data/ChIA-PET/top100_1.fq --fq2 /Users/yguo/Desktop/PPG_Data/ChIA-PET/top100_2.fq --match ATCTTATCTGAC --full ACGCGATATCTTATCTGACT --len 42
	public static void main(String[] args) {
		String fq1 = Args.parseString(args, "fq", null);
//		String fq2 = Args.parseString(args, "fq2", null);
		String match = Args.parseString(args, "match", null);	// partial fwd linker for matching
		String full = Args.parseString(args, "full", null);		// full fwd linker
		int len = Args.parseInteger(args, "len", 0);	
		int mLen = match.length();
		int idx = full.indexOf(match);		// the position of match linker inside the full linker

		ArrayList<String> l1 = CommonUtils.readTextFile(fq1);
//		ArrayList<String> l2 = CommonUtils.readTextFile(fq2);
		ArrayList<String> l1n = new ArrayList<String>();
//		ArrayList<String> l2n = new ArrayList<String>();
		
		StringBuilder sb = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		int[] labels = new int[l1.size()/4];
		for (int i=0;i<l1.size();i+=4){
			String s1 = l1.get(i+1);
			int id1 = s1.indexOf(match);
			String s1n = null;
			String q1n = null;
			if (id1>=0){
				int rightStart = id1+mLen+1;
				String q1 = l1.get(i+3);
				if (id1<=idx)
					id1 =0;
				else
					id1 -= idx;
				if (len-rightStart > id1){
					s1n = s1.substring(rightStart, len);
					q1n = q1.substring(rightStart, len);
					labels[i/4] = 1;		// 1 forward linker found, seq on right
				}
				else{
					s1n = s1.substring(0, id1);
					q1n = q1.substring(0, id1);
					labels[i/4] = 2;		// 2 forward linker found, seq on left
				}
				l1.set(i+1, s1n);
				l1.set(i+3, q1n);
			}
			else
				labels[i/4] = 0;		// 0 for linker not found
		}
		
		String match_rc = SequenceUtils.reverseComplement(match);
		int idx_rc = SequenceUtils.reverseComplement(full).indexOf(match_rc);		// the position of match linker inside the full linker

		for (int n=0;n<labels.length;n++){
			if (labels[n]!=0)		// skip if fwd linker has been found
				continue;
			int i = n*4;
			String s1 = l1.get(i+1);
			int id1 = s1.indexOf(match_rc);
			String s1n = null;
			String q1n = null;
			if (id1>=0){
				int rightStart = id1+mLen+1;
				String q1 = l1.get(i+3);
				if (id1<=idx_rc)
					id1 =0;
				else
					id1 -= idx_rc;
				if (len-rightStart > id1){
					s1n = s1.substring(rightStart, len);
					q1n = q1.substring(rightStart, len);
					labels[i/4] = 3;		// 3 revcomp linker found, seq on right
				}
				else{
					s1n = s1.substring(0, id1);
					q1n = q1.substring(0, id1);
					labels[i/4] = 4;		// 4 revcomp linker found, seq on left
				}
				l1.set(i+1, s1n);
				l1.set(i+3, q1n);
			}
			// else: 0 for linker not found
		}
				
		for (int n=0;n<labels.length;n++)
			sb.append(n).append("\t").append(labels[n]).append("\n");
		for (int n=0;n<l1.size();n++)
			sb2.append(l1.get(n)).append("\n");
		CommonUtils.writeFile(fq1+".txt", sb.toString());
		CommonUtils.writeFile(fq1+".valid.fq", sb2.toString());
	}

}
