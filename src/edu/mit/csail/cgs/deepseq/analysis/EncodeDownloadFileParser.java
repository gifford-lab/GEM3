package edu.mit.csail.cgs.deepseq.analysis;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;

public class EncodeDownloadFileParser {

	/**
	 * This java class is to parse the file.txt from ENCODE USCS download folder.
	 * 1. remove lines with "objStatus=revoked or objStatus=replaced"
	 * 
	 */
	public static void main(String[] args) {
		ArrayList<String> texts = CommonUtils.readTextFile(args[0]);
		ArrayList<String> names = new ArrayList<String>();
		ArrayList<HashMap<String,String>> infos = new ArrayList<HashMap<String,String>>();
		// parse data
		for (String s:texts){
			String[] f0 = s.split("\t");
			names.add(f0[0]);
			String[] fs = f0[1].split("; ");
			HashMap<String,String> info = new HashMap<String,String>();
			infos.add(info);
			for (String f:fs){
				String[] f1 = f.split("=");
				info.put(f1[0], f1[1]);
			}
		}
		// get all field names
		HashSet<String> fields = new HashSet<String>();
		for (HashMap<String,String> info:infos){
			fields.addAll(info.keySet());
		}
		ArrayList<String> sortedFields = new ArrayList<String>();
		sortedFields.addAll(fields);
		Collections.sort(sortedFields);
		
		// print data in table format
		// header
		StringBuilder sb = new StringBuilder();
		sb.append("#FileName");
		for (String f:sortedFields){
			sb.append("\t").append(f);
		}
		sb.append("\n");
		// data
//		String subtype = "table";
		String subtype = "fastq_5celltype_toLoad";
		SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
		Date lastLoadingDate = null;
		try {
			lastLoadingDate = format.parse("2012-03-03");
		} catch (ParseException e1) {
			e1.printStackTrace();
		}
		StringBuilder obsolete = new StringBuilder();
		for (int i=0;i<names.size();i++){
			HashMap<String,String> info = infos.get(i);
			
			// output fastq files only
			if (info.get("type").equals("fastq")){

				// remove obsolete data
				String objStatus = info.get("objStatus");
				if (objStatus!=null && (objStatus.contains("revoked") || objStatus.contains("replaced") || objStatus.contains("renamed"))){
					obsolete.append(names.get(i)).append("\t").append(objStatus).append("\n");
					continue;
				}
					
				// limit by dateUnrestricted
				String dateUnrestricted = info.get("dateUnrestricted");
				try {
					Date restrictedDate = format.parse(dateUnrestricted);
					if (restrictedDate.before(lastLoadingDate))
						continue;
				} catch (ParseException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				// limit by cell type H1-hESC
				String cell = info.get("cell");
				if (!(cell.equals("K562")||cell.equals("GM12878")||cell.equals("H1-hESC")||cell.equals("HepG2")||cell.equals("HeLa-S3")))
					continue;
				
				sb.append(names.get(i));			
				for (String key:sortedFields){
					sb.append("\t");
					if (info.containsKey(key))
						sb.append(info.get(key));
					else
						sb.append("--");
				}
				sb.append("\n");
			}
		}
		CommonUtils.writeFile(args[0]+"_revoked.txt", obsolete.toString());
		CommonUtils.writeFile(args[0]+"_"+subtype+".txt", sb.toString());
	}

}
