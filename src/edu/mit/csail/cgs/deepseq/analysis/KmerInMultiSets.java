package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class KmerInMultiSets {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// load sequence file
		HashMap<String, String> p2seqs = new HashMap<String, String>();
		ArrayList<String> strs = CommonUtils.readTextFile(Args.parseString(args, "seqs", null));
		String p2 = null;
        String[]f = null;
		for (String line: strs){
            if (line.startsWith(">")){
        		f = line.replaceFirst(">", "").split(" ");
        		if (f.length>0)
	            	p2 = f[0];
        	}
        	else{
        		p2seqs.put(p2, line);
        	}
		}
		
		// load k-mer file		
		ArrayList<String> kmers = new ArrayList<String>();
		strs = CommonUtils.readTextFile(Args.parseString(args, "kmers", null));
		for (String line: strs){
			f = line.split("\t");
			kmers.add(f[0]);
		}
		
		// load sets of points
		ArrayList<String> setNames = new ArrayList<String>();
		Vector<String> setTags=new Vector<String>();
		for(String s : args)
        	if(s.contains("set"))
        		if(!setTags.contains(s))
        			setTags.add(s);
		
		if(setTags.size()==0){
		    System.err.println("Error: No set of points provided.\n");
		    System.exit(1);
		}
        // each tag represents a file containing a subset of points
        for(String tag : setTags){
        	String name="";
        	name = tag.replaceFirst("--set", ""); 
        	setNames.add(name);
        }
        // data structure
        int[][] kmerHits = new int[kmers.size()][setNames.size()];
        int[] kmerSetHits = new int[setNames.size()];					// the collective hit count of all k-mers
        int[] setSizes = new int[setNames.size()];
        // count k-mers in each set
        for (int j =0;j<setNames.size();j++){
        	String s = setNames.get(j);
        	String filename = Args.parseString(args, "--set"+s, null);
        	ArrayList<String> ps = readPointStrings(filename);
        	setSizes[j] = ps.size();
        	for (String p:ps){
        		boolean isHit = false;		// if the sequence is hit by any k-mer
        		if (!p2seqs.containsKey(p)){
        			System.err.println(p+" sequence is not found.");
        			continue;
        		}
        		String seq = p2seqs.get(p);
        		for (int i=0;i<kmers.size();i++){
        			String km = kmers.get(i);
        			if (seq.contains(km)){
        				isHit = true;
        				kmerHits[i][j]++;
        				continue;			//skip the rc, avoid double counting
        			}
        			else{
        				seq = SequenceUtils.reverseComplement(seq);
        				if (seq.contains(km)){
            				isHit = true;
            				kmerHits[i][j]++;
            				continue;
            			}
        			}
        		}
    			if (isHit)
    				kmerSetHits[j]++;
        	}
        }
        
        //output
        StringBuilder sb =new StringBuilder();
        int k = kmers.get(0).length();
        // header
        sb.append("Kmer").append(CommonUtils.padding(k-4, " ")).append("\t");
        for (String s:setNames){
        	sb.append(s).append("\t");
        }
        sb.append("\n");
        // k-mer data
        for (int i=0;i<kmers.size();i++){
			String km = kmers.get(i);
			sb.append(km).append("\t");
			for (int j =0;j<setNames.size();j++){
				sb.append(kmerHits[i][j]).append("\t");
			}
	        sb.append("\n");
        }
        // set hit
        sb.append("SetHit").append(CommonUtils.padding(k-6, " ")).append("\t");
		for (int j =0;j<setNames.size();j++){
			sb.append(kmerSetHits[j]).append("\t");
		}
        sb.append("\n");
        sb.append("SetSize").append(CommonUtils.padding(k-7, " ")).append("\t");
		for (int j =0;j<setNames.size();j++){
			sb.append(setSizes[j]).append("\t");
		}
        sb.append("\n");
        System.out.println(sb.toString());
        CommonUtils.writeFile(Args.parseString(args, "out", null), sb.toString());
	}

	private static ArrayList<String> readPointStrings(String filename){
		ArrayList<String> ps = new ArrayList<String>();
		ArrayList<String> strs = CommonUtils.readTextFile(filename);
		for (String line: strs){
			if (line.length()==0 || line.startsWith("#"))
				continue;
			String[] f = line.split("\t");
			ps.add(f[0]);
		}
		return ps;
	}
}
