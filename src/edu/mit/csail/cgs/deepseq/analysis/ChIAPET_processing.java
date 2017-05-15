package edu.mit.csail.cgs.deepseq.analysis;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Scanner;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class ChIAPET_processing {
	// --fq1 /Users/yguo/Desktop/PPG_Data/ChIA-PET/top100_1.fq --fq2 /Users/yguo/Desktop/PPG_Data/ChIA-PET/top100_2.fq --match ATCTTATCTGAC --full ACGCGATATCTTATCTGACT --len 42
	public static void main(String[] args) {
		String fq = Args.parseString(args, "fq", null);
		String match = Args.parseString(args, "match", null);	// partial fwd linker for matching
		int len = Args.parseInteger(args, "len", 0);	
//		ArrayList<String> l1 = CommonUtils.readTextFile(fq1);
//		int[] labels = new int[l1.size()/4];
		String full = Args.parseString(args, "full", null);		// full fwd linker
		int fl = full.length();
		
		int hl = fl/2;	// half length
		String m[] = new String[hl+1];
		String m_rc[] = new String[hl+1];
		String m2[] = new String[hl+1];
		String m2_rc[] = new String[hl+1];
		for (int j=1;j<=hl;j++){
			m[j] = full.substring(j, fl);					// right part, match to the left start
			m_rc[j] = SequenceUtils.reverseComplement(m[j]);
			m2[j] = full.substring(0, fl-j);				// left part, match to the right end
			m2_rc[j] = SequenceUtils.reverseComplement(m2[j]);
		}
		
		int mLen = match.length();
		int idx = full.indexOf(match);		// the position of match linker inside the full linker
		String match_rc = SequenceUtils.reverseComplement(match);
		int idx_rc = SequenceUtils.reverseComplement(full).indexOf(match_rc);		// the position of match linker inside the full linker

		FileInputStream inputStream = null;
		Scanner sc = null;
		ArrayList<String> l1 = new ArrayList<String>();		// only hold 4 lines
		StringBuilder sb = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		CommonUtils.writeFile(fq+".txt", sb.toString());	// clean up old file
		CommonUtils.writeFile(fq+".valid.fq", sb2.toString());	
		try {
		    inputStream = new FileInputStream(fq);
		    sc = new Scanner(inputStream, "UTF-8");
		    int n = 1;
		    while (sc.hasNextLine()) {
		    	l1.clear();
		        l1.add(sc.nextLine());
		        l1.add(sc.nextLine());
		        l1.add(sc.nextLine());
		        l1.add(sc.nextLine());	// read 4 lines a time
		        
				int label = 0;
				String s1 = l1.get(1);
				int id1 = s1.indexOf(match);
				String s1n = null;
				String q1n = null;
				if (id1>=0){
					int rightStart = id1+mLen+1;
					String q1 = l1.get(3);
					if (id1<=idx)
						id1 =0;
					else
						id1 -= idx;
					if (len-rightStart > id1){
						s1n = s1.substring(rightStart, len);
						q1n = q1.substring(rightStart, len);
						label = 1;		// 1 forward linker found, seq on right
					}
					else{
						s1n = s1.substring(0, id1);
						q1n = q1.substring(0, id1);
						label = 2;		// 2 forward linker found, seq on left
					}
					if (label!=0){
						l1.set(1, s1n);
						l1.set(3, q1n);
					}
				}
				else{    // 0 for linker not found
					label = 0;		
					id1 = s1.indexOf(match_rc);
					s1n = null;
					q1n = null;
					if (id1>=0){
						int rightStart = id1+mLen+1;
						String q1 = l1.get(3);
						if (id1<=idx_rc)
							id1 =0;
						else
							id1 -= idx_rc;
						if (len-rightStart > id1){
							s1n = s1.substring(rightStart, len);
							q1n = q1.substring(rightStart, len);
							label = 3;		// 3 revcomp linker found, seq on right
						}
						else{
							s1n = s1.substring(0, id1);
							q1n = q1.substring(0, id1);
							label = 4;		// 4 revcomp linker found, seq on left
						}
						if (label!=0){
							l1.set(1, s1n);
							l1.set(3, q1n);
						}
					}
				}
				if (label<=0){
					// progressively trim to match partial linker at the edges
					for (int j=1;j<=hl;j++){
						int ml = m[j].length();

						if (s1.startsWith(m[j])){		// match
							l1.set(1, s1.substring(ml, len));
							l1.set(3, l1.get(3).substring(ml, len));
							label = 1;		// 1 forward linker found, seq on right
							break;
						}
						else if (s1.endsWith(m_rc[j])){		// match
							l1.set(1, s1.substring(0, len-ml));
							l1.set(3, l1.get(3).substring(0, len-ml));
							label = 4;		// 4 revcomp linker found, seq on left
							break;
						}
						else if (s1.startsWith(m2_rc[j])){		// match
							l1.set(1, s1.substring(ml, len));
							l1.set(3, l1.get(3).substring(ml, len));
							label = 3;		// 3 : revcomp linker found, seq on right
							break;
						}
						else if (s1.endsWith(m2[j])){		// match
							l1.set(1, s1.substring(0, len-ml));
							l1.set(3, l1.get(3).substring(0, len-ml));
							label = 2;		// 2 : forward linker found, seq on left
							break;
						}
					}	// progressive trimming
				}	// search for partial match
				
				// output
				sb.append(n).append("\t").append(label).append("\n");
				n++;
				if (sb.length()>1000000){
					CommonUtils.appendFile(fq+".txt", sb.toString());
					sb = new StringBuilder();
				}
				for (int i=0;i<4;i++)
					sb2.append(l1.get(i)).append("\n");
				if (sb2.length()>1000000){
					CommonUtils.appendFile(fq+".valid.fq", sb2.toString());
					sb2 = new StringBuilder();
				}
		    }
		    // note that Scanner suppresses exceptions
		    if (sc.ioException() != null) {
		        throw sc.ioException();
		    }
		} 
		catch (Exception e){
			e.printStackTrace();
		}
		finally {
			try{
			    if (inputStream != null) {
			        inputStream.close();
			    }
			    if (sc != null) {
			        sc.close();
			    }
			}
			catch (Exception e){
				e.printStackTrace();
			}
		}
		if (sb.length()!=0)
			CommonUtils.appendFile(fq+".txt", sb.toString());
		if (sb2.length()!=0)
			CommonUtils.appendFile(fq+".valid.fq", sb2.toString());	
	}

}
