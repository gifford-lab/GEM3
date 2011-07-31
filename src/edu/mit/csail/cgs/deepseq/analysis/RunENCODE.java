package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

public class RunENCODE {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		HashSet<String> labs = new HashSet<String>(){{
		    add("Bernstein");
		    add("Crawford");
		    add("Myers");
		    add("Snyder");
		    add("Stam");
		}};
		HashSet<String> TFtoRun = new HashSet<String>(){{
		    add("c-Myc");
//		    add("GR");
		}};
		TFtoRun.add("Input");
		
		ArrayList<String> expts = new ArrayList<String>();
		TreeSet<String> tfs = new TreeSet<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			
	        String line;
	        bin.readLine();	// skip header
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            String[] f = line.split("\t");
	            String expt = f[0];
	            String lab = f[1];
	            String cond = f[2];
	            String tf = f[3];
	            String rep = f[4];
	            String cell = f[5];
	            String genome = f[8];
	            if (genome.equalsIgnoreCase("hg19") && labs.contains(lab)  && TFtoRun.contains(tf)
	            		&& !((tf.startsWith("H3")||tf.startsWith("H4")
	            				||tf.startsWith("Pol2")
	            				||tf.startsWith("CTCF")
	            				||tf.contains("Seq")||tf.contains("seq")||tf.contains("Dnase")
	            				||tf.equals("Large-Fragment")||tf.equals("MNase")||tf.equals("Naked-DNA")
	            				||tf.equals("Mouse-IgG")||tf.equals("Control")))
	            					){
	            	expts.add(expt);
	            	tfs.add(tf);
	            }
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[0]);
            e.printStackTrace(System.err);
        }
        Collections.sort(expts);
        TreeMap<String, ArrayList<String>>exptRep = new TreeMap<String, ArrayList<String>>();
        for (String s:expts){
//        	System.out.println(s);
        	String[] f = s.split(";");
        	String name = f[0];
        	String rep = f[1];
        	if (!exptRep.containsKey(name)){
        		exptRep.put(name, new ArrayList<String>());
        	}
        	exptRep.get(name).add(s);
        }
        
        StringBuilder list_sb = new StringBuilder();
        StringBuilder gem_sb = new StringBuilder();
        for (String s:exptRep.keySet()){
        	String[] f = s.split(" ");
        	String tf = f[2];
        	if (tf.equalsIgnoreCase("Input"))
        		continue;
        	String input = s.replace(tf, "Input");
        	if (exptRep.containsKey(input)){
        		String tf_name = s.replace(' ', '_');
        		String input_name = input.replace(' ', '_');
        		list_sb.append(tf_name+"\t");
        		for (String rep:exptRep.get(s)){
        			System.out.print(rep+"\t");
        			list_sb.append(rep+"\t");
        		}
        		list_sb.deleteCharAt(list_sb.length()-1);		// remove TAB
        		list_sb.append("\n");
        		list_sb.append(input_name+"\t");
        		for (String rep:exptRep.get(input)){
        			System.out.print(rep+"\t");
        			list_sb.append(rep+"\t");
        		}
        		list_sb.deleteCharAt(list_sb.length()-1);		// remove TAB
        		list_sb.append("\n");
        		String[] ff = tf_name.split("_");
        		String run_name = ff[2]+"_"+ff[1]+"_"+ff[0];
        		gem_sb.append(run_name+"\t").append(tf_name+"\t").append(input_name+"\t--k_min 7 --k_max 13\n");
        	}
        	else{
        		for (String n:exptRep.get(s))
        			System.out.print(n+"\t");
        		System.out.print("NO INPUT");
        	}
        	System.out.println();
        }
        System.out.println("*****************************");
        System.out.println(list_sb.toString());
        System.out.println(gem_sb.toString());
	}

}
