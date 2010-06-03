package edu.mit.csail.cgs.deepseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;

/**
 * A BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * The probability is directional (i.e. stranded). 
 * Given a signed distance from the event, the BindingModel should return a relative probability of seeing a read at that distance.
 * This version of BindingModel allows support for course-grained resolutions. 
 * 
 * Not yet tested!
 * 
 * @author shaunmahony
 *
 */
public class ResolutionBindingModel {
	protected int min, max;
	protected double[] data;
	protected double[] probs;
	protected int resolution=1; //basepair resolution of the model
	protected List<Pair<Integer, Double>> bindingDistrib;

	public ResolutionBindingModel(File f){this(f,1);}
	public ResolutionBindingModel(File f, int res){
		min=0; max=0;
		resolution= res <1 ? 1 : res;
		try {
			List<Pair<Integer, Double>> bd = new LinkedList<Pair<Integer,Double>>(); 
			BufferedReader reader = new BufferedReader(new FileReader(f));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
	            	Pair<Integer,Double> p = new Pair<Integer,Double>(new Integer(words[0]), new Double(words[1]));
	            	bd.add(p);
	            }
	        }
	        bindingDistrib=bd;
	        loadData(bd);
			makeProbabilities();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public ResolutionBindingModel(List<Pair<Integer, Double>> bindingDist){this(bindingDist,1);}
	public ResolutionBindingModel(List<Pair<Integer, Double>> bindingDist, int res){
		min=0; max=0;
		bindingDistrib=bindingDist;
		loadData(bindingDist);
		makeProbabilities();
	}
	
	//Load data
	private void loadData(List<Pair<Integer, Double>> bindingDist){
		//Assumes the list is sorted//
		
		//Find max & min values first
		for(Pair<Integer, Double> p : bindingDist){
			if(p.car()<min)
				min=p.car();
			if(p.car()>max)
				max=p.car();
		}
		//Initialize arrays
		int wins = ((max-min)/resolution);
		data = new double[wins+1];
		probs = new double[wins+1];
		for(int i=0; i<=wins; i++){
			data[i]=0; probs[i]=0;
		}
		
		//Populate the data array (assumes sorted)
		int last=min-1;
		for(Pair<Integer, Double> p : bindingDist){
			int index = p.car();
			double val = p.cdr();
			//if list is not properly sorted (need to make this into an exception)
			if(index-last<0){
				System.err.println("Incorrectly sorted binding model data!"); System.exit(1);
			}
			//if unevenly spaced, smooth between values
			if(index-last>1){
				double lastVal=dataVal(last);
				double step = (val-lastVal)/(double)(index-last);
				for(int i=1; i<(index-last); i+=resolution){
					data[((last+i)-min)/resolution]=lastVal+(step*(double)i);
				}
			}
			data[(index-min)/resolution]=val;
			
			last = p.car();
		}
	}
	
	//Set a probability landscape according to the data. 
	public void makeProbabilities(){
		double totalVal=0;
		for(int i=min; i<=max; i+=resolution){
			totalVal+=dataVal(i);
		}
		for(int i=min; i<=max; i+=resolution){
			probs[(i-min)/resolution] = dataVal(i)/totalVal; 
		}
	}
	//Look up the probability corresponding to a distance
	public double probability(int distance){
		if(distance<min || distance>max){
			return(0.0);
		}else{
			return(probs[(distance-min)/resolution]);
		}
	}
	//Look up the data corresponding to a distance
	public double dataVal(int distance){
		if(distance<min || distance>max){
			return(0.0);
		}else{
			return(data[(distance-min)/resolution]);
		}
	}
	//Return the distance from the maximum value to the zero point
	public int maxShift(){
		int shift = 0; 
		double maxVal=0;
		for(int i=0; i<max/resolution; i++){
			if(data[i]>maxVal){
				shift=i; maxVal=data[i];
			}
		}return(shift*resolution);
	}
	
	//Print probs to a file
	public void printToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int i=min; i<=max; i+=resolution){
				int j = i+resolution/2;
				fout.write(j+"\t"+probability(i)+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	//Command-line interface to load a BindingModel from a file
	public static void main(String[] args){
		if(Args.parseArgs(args).contains("in")){
			String infile = Args.parseString(args, "in", null);
			String outfile = Args.parseString(args, "out", "out.model");
			File pFile = new File(infile);
			if(!pFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        
	        //File loaded, make a BindingModel
	        ResolutionBindingModel model = new ResolutionBindingModel(pFile, 1);
	        model.printToFile(outfile);
		}else{
			System.out.println("Usage: BindingModel --in filename --out outfile");
		}
	}
}
