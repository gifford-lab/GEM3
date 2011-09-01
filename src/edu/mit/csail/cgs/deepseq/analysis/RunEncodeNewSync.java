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
	
	RunEncodeNewSync(String fileName){
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
//		    add("TAF1");
//		    add("TCF12");
//		    add("TCF4");
//		    add("TFIIIC-110");
//		    add("USF-1");
//		    add("YY1");
//		    add("ZBTB33");
//		    add("c-Fos");
//		    add("c-Jun");
//		    add("p300");
			add("CTCF");
		}};
		TFtoRun.add("Input");
		
		ArrayList<String> expts = new ArrayList<String>();
		TreeSet<String> tfs = new TreeSet<String>();
		
		HashMap<String, Expt> id2expt = new HashMap<String, Expt>();
		HashMap<String, Expt> tag2input = new HashMap<String, Expt>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
			
	        String line;
	        bin.readLine();	// skip header
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            String[] f = line.split("\t");
	            String loaded = f[0];
	            if (loaded.equals("NOTLOADED"))
	            	continue;
	            String id_rep = f[1];
	            String submit = f[2];
	            String embargo = f[3];
	            String lab = f[4];
	            String group = f[5];
	            String type = f[6];
	            String cond = f[7];
	            String expt = f[8];
	            String tf = expt.replaceFirst("-ChipSeq", "");
	            if (tf.startsWith("eGFP")){
	            	String tf_f[] = expt.split("-");
	            	tf=tf_f[1]+"-"+tf_f[0];
	            }
	            String rep = f[9];
	            String cell = f[10];
	            String control = f[11];
	            String ei = f[12];
	            String url = f[13];
	            String file = f[14];
	            String size = f[15];
	            String id = String.format("ENCh-%s-%s %s %s %s", lab, group, cond, expt, cell);
	            if (ei.equals("exp")){
	            	if (tf.startsWith("H2")||tf.startsWith("H3")||tf.startsWith("H4")||tf.startsWith("Control")||tf.startsWith("Pol2"))
	            		continue;
	            	if (id2expt.containsKey(id))
	            		id2expt.get(id).addRep(rep);
	            	else{
		            	Expt e = new Expt(id, control, rep, String.format("%s_%s_%s", tf, group, cond));
	            		id2expt.put(id, e);
	            	}
	            }
	            else if (ei.equals("input")){
	            	if (tag2input.containsKey(control))
	            		tag2input.get(control).addRep(rep);
	            	else{
		            	Expt e = new Expt(id, control, rep, String.format("%s_%s_%s", tf, group, cond));
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
        exptList.addAll(id2expt.values());
        Collections.sort(exptList);
        for (Expt e:exptList){
        	expt_sb.append(String.format("%s", e.name));
        	for (String rep:e.reps)
        		expt_sb.append(String.format("\t%s;%s;bowtie_unique", e.id, rep));
        	expt_sb.append("\n");
        	if (tag2input.containsKey(e.ctrlTag))
        		gem_sb.append(String.format("%s\t%s\t%s\t--k_min 7 --k_max 13\n", e.name, e.name, tag2input.get(e.ctrlTag).name));
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
