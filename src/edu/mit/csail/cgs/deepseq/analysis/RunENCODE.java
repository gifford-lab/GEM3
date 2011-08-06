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
//		    add("c-Myc");
//		    add("GR");
//			// 0 discoveries
//		    add("BRF1");
//		    add("BRF2");
//		    add("NELFe");
//		    add("XRCC4");
//		    add("SUZ12");
//		    add("ZNF274");
//		    add("ZZZ3");
//			
////		    
////		    // 0 known motif
//		    add("BCL11A");
//		    add("BCL3");
//		    add("BAF155");
//		    add("BATF");
//		    add("BCL");
//		    add("BDP1");
//		    add("CCNT2");
//		    add("CHD2");
//		    add("CTCFL");
//		    add("HDAC2");
//		    add("HEY1");
//		    add("HMGN3");
//		    add("KAP1");
//		    add("PU.1");
//		    add("Rad21");
//		    add("SETDB1");
//		    add("SIRT6");
//		    add("SMC3");
//		    add("SP2");
//		    add("Sin3Ak-20");
//		    add("THAP1");
//		    add("TR4");
//		    add("ZNF263");
		    
//		    // else
//		    add("AP-2alpha");
//		    add("AP-2gamma");
//		    add("ATF3");
//		    add("BAF170");
//		    add("BHLHE40");
//		    add("Brg1");
//		    add("CEBPB");
//		    add("E2F4");
//		    add("E2F6");
//		    add("EBF");
//		    add("ERRA");
//		    add("Egr-1");
//		    add("FOSL2");
//		    add("FOXP2");
//		    add("GABP");
//		    add("GATA-1");
//		    add("GATA-2");
//		    add("GRp20");
//		    add("GTF2B");
//		    add("HNF4A");
//		    add("HSF1");
//		    add("IRF4");
//		    add("Ini1");
//		    add("JunD");
//		    add("Max");
//		    add("NF-E2");
//		    add("NF-YA");
//		    add("NF-YB");
//		    add("NFKB");
//		    add("NRSF");
//		    add("Nrf1");
//		    add("PAX5-C20");
//		    add("PAX5-N19");
//		    add("PGC1A");
//		    add("POU2F2");
//		    add("Pbx3");
//		    add("Pol3");
//		    add("RPC155");
//		    add("RXRA");
//		    add("SIX5");
//		    add("SP1");
//		    add("SREBP1");
//		    add("SREBP1A");
//		    add("SREBP2");
//		    add("SRF");
//		    add("STAT1");
//		    add("STAT2");
		    add("TAF1");
		    add("TCF12");
		    add("TCF4");
		    add("TFIIIC-110");
		    add("USF-1");
		    add("YY1");
		    add("ZBTB33");
		    add("c-Fos");
		    add("c-Jun");
		    add("p300");
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
	            if (genome.equalsIgnoreCase("hg19") && labs.contains(lab) && TFtoRun.contains(tf)
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
        for (String t:tfs)
        	System.out.println("add(\""+t+"\");");
        
        StringBuilder tf_sb = new StringBuilder();
        StringBuilder input_sb = new StringBuilder();
        StringBuilder gem_sb = new StringBuilder();
        for (String s:exptRep.keySet()){
        	String[] f = s.split(" ");
        	String tf = f[2];
        	if (tf.equalsIgnoreCase("Input"))
        		continue;
        	String input = s.replace(tf, "Input");
        	if (exptRep.containsKey(input)){
        		String tf_name = s.replace('_', '-');
        		tf_name = tf_name.replace(' ', '_');
        		String input_name = input.replace('_', '-');
        		input_name = input_name.replace(' ', '_');
        		tf_sb.append(tf_name+"\t");
        		for (String rep:exptRep.get(s)){
        			System.out.print(rep+"\t");
        			tf_sb.append(rep+"\t");
        		}
        		tf_sb.deleteCharAt(tf_sb.length()-1);		// remove TAB
        		tf_sb.append("\n");
        		input_sb.append(input_name+"\t");
        		for (String rep:exptRep.get(input)){
        			System.out.print(rep+"\t");
        			input_sb.append(rep+"\t");
        		}
        		input_sb.deleteCharAt(input_sb.length()-1);		// remove TAB
        		input_sb.append("\n");
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
        System.out.println(tf_sb.toString());
        System.out.println(input_sb.toString());
        System.out.println(gem_sb.toString());
	}

}
