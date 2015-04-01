package edu.mit.csail.cgs.deepseq.discovery.kmer;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public class KmerSet{
	private ArrayList<Kmer> kmers;
	public int posSeqCount=-1;
	public int negSeqCount=-1;
	public double ksmThreshold=0;
	private TreeSet<Integer> clusterIds;
	
	public KmerSet (File file){
		kmers = new ArrayList<Kmer>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	        String line = bin.readLine().trim();					// first line, pos_seq# / neg_seq#
	        line = line.substring(1,line.length());			//remove # sign
	        String[] f = line.split("/");
	        posSeqCount = Integer.parseInt(f[0]);
	        negSeqCount = Integer.parseInt(f[1]);
	        
	        line = bin.readLine().trim();							// 2nd line, KCM threshold
	        line = line.substring(1,line.length());			//remove # sign
	        ksmThreshold = Double.parseDouble(line);

	        clusterIds = new TreeSet<Integer>();
	        while((line = bin.readLine()) != null) { 
	        	if (line.startsWith("#"))
	        		continue;
	            line = line.trim();
	            Kmer kmer = Kmer.fromString(line);
	            kmers.add(kmer);
	            clusterIds.add(kmer.getClusterId());
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+file.getName());
            e.printStackTrace(System.err);
        }
        kmers.trimToSize();
	}

	public ArrayList<Kmer> getKmers (int clusterId){
		ArrayList<Kmer> selected = new ArrayList<Kmer>();
		for (Kmer km: kmers){
			if (km.getClusterId()==clusterId)
				selected.add(km);
		}
		return selected;
	}
	public Set<Integer> getClusterIds (){
		return clusterIds;
	}
	
	public static void main0(String[] args){
		KmerSet kc = new KmerSet(new File(args[0]));
		return;
	}
}
