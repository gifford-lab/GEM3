package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;

public class EncodeNameMapping {

	public static void main(String[] args) {
	// Command line: C:\Data\ENCODE\MotifCompare\name_map.txt C:\Data\ENCODE\MotifCompare\top_encode_PFM.txt
//		mapToFactorGroup(args);
		
	// Command line: C:\Data\ENCODE\MotifCompare\info\ENCODE_name_mapping_Pouya.txt C:\Data\ENCODE\MotifCompare\Aug_runs\top_encode_PFM_toPouya.txt
//		mapToExpName(args);
	// Command line: C:\Data\ENCODE\MotifCompare\info\ENCODE_name_mapping_Pouya.txt C:\Data\ENCODE\MotifCompare\encode_PFM_1.txt
		mapToExpNameNewSync(args);
	}
	public static void mapToFactorGroup(String[] args) {
		// TODO Auto-generated method stub
		HashMap<String, String> map = new HashMap<String, String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			
	        String line;
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	String[] f = line.split("\t");
	        	map.put(f[1], f[0]);
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[0]);
            e.printStackTrace(System.err);
        }
        TreeSet<String> newTF = new TreeSet<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[1])));
			
	        String line;
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	if (!line.startsWith("DE "))
	        		continue;
	        	String[] f = line.substring(3).split("_");
	        	String tf = f[0];
	        	if (!map.containsKey(tf))
	        		newTF.add(tf);
	        	else
	        		System.out.println(map.get(tf)+"\t"+line.substring(3));
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[1]);
            e.printStackTrace(System.err);
        }
        System.out.println("**************************");
        for (String tf:newTF){
        	System.out.println("\t"+tf);
        }
	}
	
	
	// Command line: C:\Data\ENCODE\MotifCompare\ENCODE_name_mapping_Pouya.txt C:\Data\ENCODE\MotifCompare\top_encode_PFM_toPouya.txt
	public static void mapToExpName(String[] args) {
		// TODO Auto-generated method stub
		HashMap<String, String> map = new HashMap<String, String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			
	        String line;
	        bin.readLine();
	        bin.readLine();
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	String[] f = line.split("\t");
	        	String name = f[1].trim();
	        	String[] t = name.split("_");
	        	String pi = t[2].substring(7);
	        	String key;
	        	if (t[0].equals("AP-2alpha"))				// debugging line
	        		pi=pi;
	        	if (t.length==5 || (t.length==6 && (t[5].startsWith("IgG")||t[5].startsWith("v041610")||t[5].startsWith("C-8")
	        			||t[5].startsWith("PCR2x")||t[5].startsWith("M-17")||t[5].startsWith("pravastatin"))))
	        		key= t[0]+"_"+t[1]+"_"+pi;
	        	else
	        		key= t[0]+"_"+t[1]+"-"+t[5]+"_"+pi;
	        	if (map.containsKey(key))
	        		System.out.println("Duplicate Key: "+key+" => "+map.get(key) + "\t" + name);
	        	else
	        		map.put(key, name);
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[0]);
            e.printStackTrace(System.err);
        }
        TreeSet<String> newTF = new TreeSet<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[1])));
			
	        String line;
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	if (!line.startsWith("DE "))
	        		continue;
	        	String name = line.substring(3);
	        	String[] t = name.split("_");
	        	String tf = t[0];
	        	if (tf.equals("GATA-1")) tf="GATA1";
	        	if (tf.equals("GATA-2")) tf="GATA2";
	        	if (tf.equals("c-Myc")) tf="Myc";
	        	if (tf.equals("c-Jun")) tf="Jun";
	        	if (tf.equals("c-Fos")) tf="Fos";
	        	if (tf.equals("SREBP1A")) tf="SREBP1";
	        	if (tf.equals("USF-1")) tf="USF1";
	        	String key = tf+"_"+t[1]+"_"+t[2];
	        	
	        	if (!map.containsKey(key))
	        		newTF.add(name+"\t"+key);
	        	else
	        		System.out.println(name+"\t"+map.get(key));
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[1]);
            e.printStackTrace(System.err);
        }
        System.out.println("\n******* These need to be fixed manually ********");
        for (String tf:newTF){
        	System.out.println(tf);
        }
	}
	

	// Command line: C:\Data\ENCODE\MotifCompare\info\ENCODE_name_mapping_Pouya.txt C:\Data\ENCODE\MotifCompare\encode_PFM_1.txt
	public static void mapToExpNameNewSync(String[] args) {
		// TODO Auto-generated method stub
		HashMap<String, String> map = new HashMap<String, String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
			
	        String line;
	        bin.readLine();
	        bin.readLine();
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	String[] f = line.split("\t");
	        	String name = f[1].trim();
	        	String[] t = name.split("_");
	        	String pi = t[2].substring(7);
	        	String key;
	        	if (t[0].equals("EBF"))				// debugging line
	        		pi=pi;
	        	if (t.length==6 && t[5].contains("-v041610"))
	        		t[5] = t[5].substring(0, t[5].indexOf("-v041610"));
	        	if (t.length==5 || (t.length==6 && (t[5].startsWith("IgG")||t[5].startsWith("v041610")||t[5].startsWith("C-8")
	        			||t[5].startsWith("PCR2x")||t[5].startsWith("M-17")||t[5].startsWith("pravastatin")||t[5].startsWith("M33")
	        			||t[5].startsWith("estrogen"))))
	        		key= t[0]+"_"+t[1]+"_"+pi;
	        	else
	        		key= t[0]+"_"+t[1]+"-"+t[5]+"_"+pi;
//	        	if (!key.equals("ATF3_K562_Myers") && t[0].equals("BCLAF1"))				// debugging line       	
	        	if (t[0].equals("EBF"))				// debugging line
	        		pi=pi;
	        	
	        	if (map.containsKey(key))
	        		System.out.println("Duplicate Key: "+key+" => "+map.get(key) + "\t" + name);
	        	else
	        		map.put(key, name);
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[0]);
            e.printStackTrace(System.err);
        }
        
        // Some name that are not consistent by splitting into cell/condition/etc, directly map
        HashMap directMap = new HashMap<String, String>(){{
        	put("ATF3-v041610.1_Myers_K562", "ATF3_K562_encode-Snyder_seq_hsa");
        	put("c-Myc_Crawford_MCF-7-estrogen", "Myc_MCF-7_encode-Crawford_seq_hsa_estrogen");
        	put("CTCF_Crawford_MCF-7-estrogen", "CTCF_MCF-7_encode-Crawford_seq_hsa_estrogen");
        	put("EBF-PCR1x_Myers_GM12878", "EBF1_GM12878_encode-Myers_seq_hsa_C-8");
        	put("FOS-eGFP_White_K562", "Fos_K562_encode-White_seq_hsa_eGFP");
        	put("FOSL1--SC-183-v041610.1_Myers_K562", "FosL1_K562_encode-Myers_seq_hsa_v041610.1-SC-183");
        	put("FOSL2-v041610.1_Myers_HepG2", "FosL2_HepG2_encode-Myers_seq_hsa_v041610.1");
        	put("GR-PCR2x_Myers_A549-DEX-100nM", "GR_A549_encode-Myers_seq_hsa_DEX-100nM-PCR2x");
        	put("JunB-eGFP_White_K562", "JunB_K562_encode-White_seq_hsa_eGFP");
        	put("JunD-eGFP_White_K562", "JunD_K562_encode-White_seq_hsa_eGFP");
        	put("NR4A1-eGFP_White_K562", "NR4A1_K562_encode-White_seq_hsa_eGFP");
        	put("PAX5-C20-PCR1x_Myers_GM12878", "PAX5_GM12878_encode-Myers_seq_hsa_C20");
        	put("PAX5-C20-v041610.1_Myers_GM12891", "PAX5_GM12891_encode-Myers_seq_hsa_v041610.1-C20");
        	put("PAX5-C20-v041610.1_Myers_GM12892", "PAX5_GM12892_encode-Myers_seq_hsa_v041610.1-C20");
        	put("PAX5-N19-PCR1x_Myers_GM12878", "PAX5_GM12878_encode-Myers_seq_hsa_N19");
        	put("Egr-1-PCR2x_Myers_GM12878", "Egr-1_GM12878_encode-Myers_seq_hsa_PCR2x");
        	put("Egr-1-v041610.1_Myers_GM12878", "Egr-1_GM12878_encode-Myers_seq_hsa_v041610.1");
        	put("Egr-1-v041610.2_Myers_H1-hESC", "Egr-1_H1-hESC_encode-Myers_seq_hsa_v041610.2");
        	put("GATA2-eGFP_White_K562", "GATA2_K562_encode-White_seq_hsa_eGFP");
        	put("GATA2--SC-267-PCR1x_Myers_K562", "GATA2_K562_encode-Myers_seq_hsa_CG2-96");
        	put("FOXA1--SC-101058-v041610.1_Myers_HepG2", "FOXA1_HepG2_encode-Myers_seq_hsa_v041610.1-SC-101058");
        	put("HDAC8-eGFP_White_K562", "HDAC8_K562_encode-White_seq_hsa_eGFP");
        	put("NRSF-v041610.2_Myers_H1-hESC", "NRSF_H1-hESC_encode-Myers_seq_hsa_v041610.2");
        	put("p300-v041610.1_Myers_HepG2", "p300_HepG2_encode-Myers_seq_hsa_v041610.1");
        	put("SRF-v041610.1_Myers_GM12878", "SRF_GM12878_encode-Myers_seq_hsa_v041610.1");
        	put("TAF1-v041610.2_Myers_GM12892", "TAF1_GM12892_encode-Myers_seq_hsa_v041610.2");
        	put("TAF1-v041610.2_Myers_H1-hESC", "TAF1_H1-hESC_encode-Myers_seq_hsa_v041610.2");
        	put("TAF1-v041610.1_Myers_K562", "TAF1_K562_encode-Myers_seq_hsa_v041610.1");
        	put("YY1-v041610.2_Myers_K562", "YY1_K562_encode-Myers_seq_hsa_v041610.2");
        	put("ZBTB33-v041610.1_Myers_HepG2", "ZBTB33_HepG2_encode-Myers_seq_hsa_v041610.1");
        }};
        	
        // ENCODE experiment not found in Pouya's list, have been checked manually
		HashSet<String> notFound_checked = new HashSet<String>(){{
		    add("FOSL1_H1-hESC_Myers");
		    add("FOXP2_SK-N-MC_Myers");
		    add("FOXP2_PFSK-1_Myers");
		    add("GRp20_HepG2-forskolin_Snyder");
		    add("MEF2A_K562_Myers");
		    add("NRSF_SK-N-SH_Myers");
		    add("Pol3_GM12878_Snyder");
		    add("Pol3_K562_Snyder");
		    add("YY1_HepG2_Myers");
		    add("Jun_GM12878_Snyder");
		    add("YY1_HCT-116_Myers");
		    add("ZBTB33_HCT-116_Myers");
		    add("Myc_GM12878_Snyder");
		}};
        TreeSet<String> newTF = new TreeSet<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(args[1])));
			
	        String line;
	        while((line = bin.readLine()) != null) { 
	        	line = line.trim();
	        	if (!line.startsWith("DE "))
	        		continue;
	        	String name = line.substring(3);
	        	String[] t = name.split("_");
	        	String tf = t[0];
	        	String directKey = tf+"_"+t[1]+"_"+t[2];
	        	if (directMap.containsKey(directKey)){
	        		System.out.println(name+"\t"+directMap.get(directKey));
	        		continue;
	        	}
	        	
	        	if (tf.contains("-v041610"))
	        		tf = tf.substring(0, tf.indexOf("-v041610"));
	        	if (tf.contains("--"))
	        		tf = tf.substring(0, tf.indexOf("--"));
	        	tf = tf.replace("-PCR1x", "");
	        	tf = tf.replace("-PCR2x", "");
	        	
	        	if (tf.equals("GATA-1")) tf="GATA1";
	        	if (tf.equals("GATA-2")) tf="GATA2";
	        	if (tf.equals("c-Myc")) tf="Myc";
	        	if (tf.equals("c-Jun")) tf="Jun";
	        	if (tf.equals("c-Fos")) tf="Fos";
	        	if (tf.equals("SREBP1A")) tf="SREBP1";
	        	if (tf.equals("USF-1")) tf="USF1";
	        	if (tf.equals("EBF") && t[1].equals("Myers")) tf="EBF";
	        	
	        	t[2] = t[2].replace("-pravastatin", "");
	        	
	        	String key = tf+"_"+t[2]+"_"+t[1];
	        	
	        	if (!map.containsKey(key)){
	        		if (!notFound_checked.contains(key))
	        			newTF.add(name+"\t"+key);
	        	}
	        	else
	        		System.out.println(name+"\t"+map.get(key));
	        }
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+args[1]);
            e.printStackTrace(System.err);
        }
        System.out.println("\n******* These need to be fixed manually ********"+newTF.size());
        for (String tf:newTF){
        	System.out.println(tf);
        }
	}
}
