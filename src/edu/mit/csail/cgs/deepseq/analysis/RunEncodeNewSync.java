package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;

public class RunEncodeNewSync {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		RunEncodeNewSync r = new RunEncodeNewSync(args[0]);
	}
	
	/**
	 * Take Shaun's ENCODE (new sync) download.txt file as input<br>
	 * Generate pair of ChIP-Seq IP/CTRL pairs for GEM script to process, together with file to group replicates to one single experiments.
	 * @param fileName
	 */
	RunEncodeNewSync(String fileName){
		HashMap<String, Expt> tfid2expt = new HashMap<String, Expt>();
		HashMap<String, Expt> tag2input = new HashMap<String, Expt>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
			
	        String line;
	        bin.readLine();	// skip header
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            String[] f = line.split("\t");
	            String id_rep = f[0];
	            String submit = f[1];
	            String embargo = f[2];
	            String lab = f[3];
	            String pi = f[4];
	            String type = f[5];
	            String cond = f[6];
	            String expt = f[7];
	            String tf = expt.replaceFirst("-ChipSeq", "");
	            if (tf.startsWith("eGFP")){
	            	String tf_f[] = expt.split("-");
	            	tf=tf_f[1]+"-"+tf_f[0];
	            }
	            String rep = f[8];
	            String cell = f[9];
	            if (!cell.equals("H1-hESC"))
	            	continue;
	            String control = lab+"_"+pi+"_"+f[10];	// add lab+PI name to distinguish inputs for the same cell from different labs
	            String ei = f[11];		// exp identifier: exp or input
	            String url = f[12];
	            String file = f[13];
	            String size = f[14];
	            String id = String.format("ENCh-%s-%s %s %s %s", lab, pi, cond, expt, cell);		//TF ID
	            if (ei.equals("exp")){
	            	if (tf.startsWith("H2")||tf.startsWith("H3")||tf.startsWith("H4")||tf.startsWith("Control"))
	            		continue;
	            	if (tfid2expt.containsKey(id))
	            		tfid2expt.get(id).addRep(rep);
	            	else{
		            	Expt e = new Expt(id, control, rep, String.format("%s_%s_%s", tf, pi, cond));
	            		tfid2expt.put(id, e);
	            	}
	            }
	            else if (ei.equals("input")){
	            	if (tag2input.containsKey(control))
	            		tag2input.get(control).addRep(rep);
	            	else{
		            	Expt e = new Expt(id, control, rep, String.format("%s_%s_%s", tf, pi, cond));
		            	tag2input.put(control, e);
	            	}
	            }
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+fileName);
            e.printStackTrace(System.err);
        }

        StringBuilder expt_sb = new StringBuilder();
        StringBuilder input_sb = new StringBuilder();
        StringBuilder gem_sb = new StringBuilder();
        ArrayList<Expt> exptList = new ArrayList<Expt>();
        exptList.addAll(tfid2expt.values());
        Collections.sort(exptList);
        for (Expt e:exptList){
        	expt_sb.append(String.format("%s", e.name));
        	for (String rep:e.reps)
        		expt_sb.append(String.format("\t%s;%s;bowtie_unique", e.id, rep));
        	expt_sb.append("\n");
        	if (tag2input.containsKey(e.ctrlTag))
        		gem_sb.append(String.format("%s\t%s\t%s\t--k_min 6 --k_max 13\n", e.name, e.name, tag2input.get(e.ctrlTag).name));
        	else
        		System.err.println(e.name+"\t"+e.ctrlTag+"\t No control!");
        }
        
        exptList.clear();
        exptList.addAll(tag2input.values());
        Collections.sort(exptList);
        for (Expt e:exptList){
        	input_sb.append(String.format("%s", e.name));
        	for (String rep:e.reps)
        		input_sb.append(String.format("\t%s;%s;bowtie_unique", e.id, rep));
        	input_sb.append("\n");
        }
        
//        gem_sb.append(run_name+"\t").append(tf_name+"\t").append(input_name+"\t--k_min 7 --k_max 13\n");

        
        System.out.println("*****************************");
        System.out.println(expt_sb.toString());
        System.out.println("*****************************");
        System.out.println(input_sb.toString());
        System.out.println("*****************************");
        System.out.println(gem_sb.toString());
	}
	class Expt implements Comparable<Expt>{
		String id;
		ArrayList<String> reps = new ArrayList<String>();
		String ctrlTag;
		String name;
		Expt(String id, String ctrlTag, String rep,  String name){
			this.id = id;
			this.ctrlTag = ctrlTag;
			this.name = name;
			reps.add(rep);
		}
		void addRep(String rep){
			reps.add(rep);
		}
		public int compareTo(Expt e) {
			return name.compareTo(e.name); // descending
		}
	}
}
