package edu.mit.csail.cgs.deepseq.analysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class ConsolidateEncodeInteractions {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ArrayList<String> texts = CommonUtils.readTextFile(args[0]);
		HashMap<String, ArrayList<Integer>> allPair2positions = new HashMap<String, ArrayList<Integer>>();
		
		HashMap<String, ArrayList<String>> entryByCell = new HashMap<String, ArrayList<String>> ();
		for (String line:texts){
			if (line.length()==0)
				continue;
			String[] f = line.split("\t");
			String cell = f[6];
			if (!entryByCell.containsKey(cell))
				entryByCell.put(cell, new ArrayList<String>());
			entryByCell.get(cell).add(line);
		}
		
		for (String cell:entryByCell.keySet()){
			HashMap<String, ArrayList<Integer>> pair2positions = new HashMap<String, ArrayList<Integer>>();
			for (String line:entryByCell.get(cell)){
				String[] f = line.split("\t");
				String tf1 = f[0];
				if (tf1.contains(":"))
					tf1=tf1.substring(0,tf1.indexOf(":"));
				String tf2 = f[1];
				if (tf2.contains(":"))
					tf2=tf2.substring(0,tf2.indexOf(":"));
				int pos = Integer.parseInt(f[3]);

				String pair;
				if (tf1.compareTo(tf2)<0)
					pair = tf1+"*"+tf2;
				else{
					pair = tf2+"*"+tf1;
//					pos = -pos;		// comment out, report distance
				}
				if (!pair2positions.containsKey(pair))
					pair2positions.put(pair, new ArrayList<Integer>());
				pair2positions.get(pair).add(pos);
			}
						
			ArrayList<String> sortedPairs = new ArrayList<String>();
			sortedPairs.addAll(pair2positions.keySet());
			if (!sortedPairs.isEmpty()){
				Collections.sort(sortedPairs);
				for (String p:sortedPairs){
					System.out.print(p.replace('*', '\t'));
					System.out.print('\t');
					ArrayList<Integer> positions = pair2positions.get(p);
					int[] elements = StatUtil.sortByOccurences(positions).car();
					System.out.print(elements[elements.length-1]+"\t");
					for (int position:positions)
						System.out.print(position+" ");
					System.out.print('\n');
				}
				System.out.println("----------------------------------");
				System.out.println(cell+"\t"+sortedPairs.size());
				System.out.println();
				
				// store to allInteraction table

				for (String p:sortedPairs){
					if (!allPair2positions.containsKey(p))
						allPair2positions.put(p, new ArrayList<Integer>());
					allPair2positions.get(p).addAll(pair2positions.get(p));
				}
			}
		}
		
		// output all cell types 
		ArrayList<String> sortedPairs = new ArrayList<String>();
		sortedPairs.addAll(allPair2positions.keySet());
		if (!sortedPairs.isEmpty()){
			Collections.sort(sortedPairs);
			for (String p:sortedPairs){
				System.out.print(p.replace('*', '\t'));
				System.out.print('\t');
				ArrayList<Integer> positions = allPair2positions.get(p);
				int[] elements = StatUtil.sortByOccurences(positions).car();
				System.out.print(elements[elements.length-1]+"\t");
				for (int position:positions)
					System.out.print(position+" ");
				System.out.print('\n');
			}
			System.out.println("----------------------------------");
			System.out.println("All cell type\t"+sortedPairs.size());
			System.out.println();
		}
	}

}
