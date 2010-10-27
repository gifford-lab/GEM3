/**
 * 
 */
package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

/**
 * @author Yuchun Guo
 * This clsss is to hold some frequently used static methods
 *
 */
public class CommonUtils {
	/*
	 * load text file in CGS Point format
	 * chr:coord, e.g. 1:234234
	 */
	public static ArrayList<Point> loadCgsPointFile(String filename, Genome genome) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				Region point = Region.fromString(genome, line);
				if (point!=null)
					points.add(new Point(genome, point.getChrom(),point.getStart()));
			}
		}
		catch(IOException ioex) {
			//logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				//logger.error("Error closing buffered reader", ioex2);
			}			
		}
		return points;
	}
	
	public static ArrayList<Region> loadRegionFile(String fname, Genome gen){
		ArrayList<Region> rset = new ArrayList<Region>();
		try{
			File rFile = new File(fname);
			if(!rFile.isFile()){
				System.err.println("\nThe region file is not found: "+fname+"!");
				System.exit(1);
			}
	        BufferedReader reader = new BufferedReader(new FileReader(rFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
            	Region r = Region.fromString(gen, words[0]);
            	if (r!=null)
            		rset.add(r);
            }
	        reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return rset;
	}
	public static String timeElapsed(long tic){
		return timeString(System.currentTimeMillis()-tic);
	}
	private static String timeString(long length){
		float sec = length/1000F;
		return sec>60?
			String.format("%.1f",sec/60)+" min":
			String.format("%.1f",sec)+" sec";
	}
	public static String getDateTime() {
        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
        Date date = new Date();
        return dateFormat.format(date);
    }
	public static void printArray(double[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static void printArray(int[]array, String msgBefore, String msgAfter){
		System.out.print(msgBefore);
		System.out.print(arrayToString(array));
		System.out.print(msgAfter);
	}
	public static String arrayToString(int[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%d\t",array[i]));
        }
        output.append(String.format("%d",array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%.2f\t",array[i]));
        }
        output.append(String.format("%.2f",array[array.length-1]));
        return output.toString();
	}
	public static String arrayToString(double[] array, int digit){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<array.length-1;i++){
        	output.append(String.format("%."+digit+"f\t",array[i]));
        }
        output.append(String.format("%."+digit+"f",array[array.length-1]));
        return output.toString();
	}	
	public static String matrixToString(double[][] matrix, int digit, String[] names){
        StringBuilder output = new StringBuilder();
        for (int i=0;i<matrix.length;i++){
        	if (names!=null)
        		output.append(names[i]).append("\t");
        	output.append(arrayToString(matrix[i], digit)).append("\n");
        }
        return output.toString();
	}
	
	public static void writeFile(String fileName, String text){
		try{
			FileWriter fw = new FileWriter(fileName, false); //new file
			fw.write(text);
			fw.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
//		System.out.println("File was written to "+fileName);
	}
    /*
     * parse motif string, version, threshold, species_id from command line
     */
    public static Pair<WeightMatrix, Double> loadPWM(String args[], int orgId){
    	Pair<WeightMatrix, Double> pair=null;
    	try {   
            String motifString = Args.parseString(args, "motif", null);
            if (motifString!=null){
      	      String motifVersion = Args.parseString(args, "version", null);
      	      double motifThreshold = Args.parseDouble(args, "motifThreshold", -1);
      	      if (motifThreshold==-1){
      	    	  System.err.println("No motif threshold was provided, default=10.0 is used.");
      	    	  motifThreshold = 10.0;
      	      }
      		  int motif_species_id = Args.parseInteger(args, "motif_species_id", -1);
      		  int wmid = WeightMatrix.getWeightMatrixID(motif_species_id!=-1?motif_species_id:orgId, motifString, motifVersion);
      		  WeightMatrix motif = WeightMatrix.getWeightMatrix(wmid);
      		  pair = new Pair<WeightMatrix, Double>(motif, motifThreshold);
            }
          } 
          catch (NotFoundException e) {
            e.printStackTrace();
          }  
  		  return pair;
    }
}
