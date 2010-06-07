package edu.mit.csail.cgs.deepseq.multicond;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.regex.Matcher;


/**
 * Class ParametersGrid creates a grid of the values of parameters
 * Parameter Names: Array of Strings representing the names of the parameters
 * Parameter Values: Array of Strings representing the values of these parameters
 * <u>Note</u>: Both parameter names and values should be in order. I.e. the first element
 * of <tt>paramNamesArg</tt> should represent the first element of <tt>paramValsArg</tt> and so on.
 * Each row of the <tt>paramValsArg</tt> array can explicitly give the values for this parameter or (in the case
 * of numeric parameters) give it in the form: <<tt> initial_value : step : final_value </tt> ><br> 
 * E.g.:
 * <pre>
 * String[] paramNamesArg = {"operation", "s",                   "class",   "lambda"};
 * String[] paramValsArg = {"*, +",       "0, 1:1:10, 11, 12",   "I, X",    "-1.5:0.1:2.5"};
 * </pre>
 * So, if for example we had:
 * <pre>
 * String[] paramNamesArg = {"x", "y"};
 * String[] paramValsArg = {"1, 2", "a, b, c"};
 * </pre>
 * The combinations would look like:
 * <pre>
 * 1	a
 * 1	b
 * 1	c
 * 2	a
 * 2	b
 * 2	c
 * </pre>
 * @author geopapa
 *
 */
public class ParametersGrid {
		
	/**
	 * E.g.: <tt>32</tt>, <tt>-12</tt>, <tt>2.3</tt>, <tt>-3.34</tt> 
	 */
	private static String numberPattString = "-?\\d*(?:\\.\\d+)?";
	private static Pattern numberPatern = Pattern.compile(numberPattString); 
	
	/**
	 * E.g.: <tt>1:1:10</tt>, <tt>1.3:.1:2.3</tt>, <tt>-2.2:0.2:3.2</tt>
	 */
	private static String stepPattString = "-?\\d*(?:\\.\\d+)?" + "\\s*:\\s*" + "-?\\d*(?:\\.\\d+)?" + "\\s*:\\s*" + "-?\\d*(?:\\.\\d+)?";
	private static Pattern stepPattern = Pattern.compile(stepPattString);
	
	/**
	 * Number of possible combinations of the different variables
	 */
	private int numCombins = 1;
	
	/**
	 * Names of the parameters
	 */
	private List<String> paramNames;
	
	/**
	 * Values of the parameters 
	 */
	private List[] paramValues;
	
	/**
	 * Possible combinations of the different variables
	 */
	private List[] combins; 

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		String[] paramNamesArg = {"operation", "s", "", "lambda"};
		
		String[] paramValsArg = {
				"*, +",
				"0, 1:.5:2, 3", //"1:1:3",   // start = 1, end = 3 in steps of 1
				"I, X",
				"-1.5 : 0.5 : 0" //"1:1:4"		// start = 1, end = 4 in steps of 1		
								};
		
						          /* { "*, +", "1, 2, 3", "I, X", "a, b, c, d"};*/
		/*
		List[] paramValsArg = new ArrayList[paramValsStr.length];
		for(int i = 0; i < paramValsArg.length; i++) {
			paramValsArg[i] = new ArrayList(Arrays.asList(paramValsStr[i]));
		}
		*/
		
		ParametersGrid grid = new ParametersGrid(paramNamesArg, paramValsArg);
		grid.evalCombinations();
		grid.printCombinations(System.out);
	}//end of main method
	
	/**
	 * 
	 * @param paramNamesArg Names of the parameters
	 * @param paramValsArg Values of the parameters
	 */
	public ParametersGrid(String[] paramNamesArg, String[] paramValsArg)
	{
		/**********************/
		// Exceptions Section // 
		/**********************/
		
		if( paramValsArg.length != paramNamesArg.length )
			throw new IllegalArgumentException("The length of the array of lists has to be equal with the sixe of the list paramNames");
		
		if(paramNamesArg.length == 0 || paramValsArg.length == 0)
			throw new IllegalArgumentException("There has to be at least one parameter examined");
		
		for( String s: paramValsArg )
			if( s.length() == 0) { throw new IllegalArgumentException("At least one of the parameters has no values"); }
		
		paramNames = Arrays.asList(paramNamesArg);
		paramValues = evalParamValues(paramValsArg); 
		
		for(List c:paramValues)
			numCombins *= c.size();
		
		
		combins = new ArrayList[numCombins];
		for(int i = 0; i < combins.length; i++) {
			combins[i] = new ArrayList(); 
		}
		
	}//end of ParametersGrid constructor
	
	
	/**
	 * Isolates the tokens of each element of <tt>paramValsArg</tt> and creates
	 * a List containing the values of the parameters in an expanded form 
	 * @param paramValsArg Values of the parameters
	 * @return
	 */
	private List[] evalParamValues(String[] paramValsArg)
	{
		List[] paramVals = new ArrayList[paramValsArg.length];
		for(int i = 0; i < paramVals.length; i++)
			paramVals[i] = new ArrayList();
			
		for(int i = 0; i < paramValsArg.length; i++)
		{
			String currParamVals = paramValsArg[i];
			
			String[] tokens = currParamVals.split("\\s*,\\s*");
			for(String tok : tokens)
			{
				Matcher m = stepPattern.matcher(tok);
				if(m.matches())
				{
					paramVals[i].addAll(unfoldNumberSequence(tok)); // method that creates
				}
				else
				{
					paramVals[i].add(tok);
				}
			}
			
		}
		
		return paramVals;
	}//end of evalParamValues method
	
	
	/**
	 * Unfolds the token when it is on the form: <<tt> initial_value : step : final_value </tt> >
	 * @param token
	 * @return
	 */
	private List<String> unfoldNumberSequence(String token)
	{
		/**********************/
		// Exceptions Section // 
		/**********************/
		String[] tokens = token.split("\\s*:\\s*");
		
		if(tokens.length > 3)
			throw new IllegalArgumentException("Short numbering format has to be in the form: start:step:end");
		
		if( (Double.parseDouble(tokens[0])-Double.parseDouble(tokens[2]))*Double.parseDouble(tokens[1]) > 0 )
			throw new IllegalArgumentException("When start < end => step has to be > 0 and when start > end => step has to be < 0 ");

		
		/******************/
		// Method's Body // 
		/*****************/
		List<String> expandedList = new ArrayList<String>(); 
		Double[] numbers = new Double[tokens.length];
		
		for(int i = 0; i < numbers.length; i++)
			numbers[i] = Double.parseDouble(tokens[i]);		
		
		double start = numbers[0];
		double end = numbers[2];
		double step = numbers[1];
		
		double currIndex = start;
		while(currIndex <= end)
		{
			expandedList.add(Double.toString(currIndex));
			currIndex += step;
		}
		
		return expandedList;
	}//end of unfoldNumberSequence method
	
	
	/**
	 * Evaluates the combinations of the values of the parameters
	 */
	public void evalCombinations()
	{		
		/******************/
		// Method's Body // 
		/*****************/
		int i, j, counter;
		int numIters = numCombins;
		
		for(i = 0; i < paramValues.length; i++)
		{			
			List c = paramValues[i];
		
			counter = 0;
			numIters /= c.size();
			j = 0;
			while( j < c.size() || counter < numCombins )
			{	
				if( j == c.size() ) { j = 0; }
				
				for(int k = 0; k < numIters; k++)
					combins[counter++].add(c.get(j));
				
				j++;
			}
		}

	}// end of evalCombinations method
	
	
	/**
	 * 
	 * @return the combinations
	 */
	public List[] getCombinations()
	{
		return combins;
	}//end of getCombinations method
	
	
	/**
	 * Prints the combinations in the specified output stream
	 * @param os
	 */
	public void printCombinations(OutputStream os)
	{
		BufferedWriter bw = null; 
		
		try {
			OutputStreamWriter osw = new OutputStreamWriter(os);
			bw = new BufferedWriter(osw);
			
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < paramNames.size()-1; i++)
				sb.append(paramNames.get(i) + "\t");
			
			sb.append(String.format("%s%n", paramNames.get(paramNames.size()-1)));
			
			bw.write(sb.toString());
			
			for(int i = 0; i < combins.length; i++)
			{
				List combin = combins[i];
				
				sb = new StringBuilder();
				for(int k = 0; k < combin.size()-1; k++)
					sb.append(combin.get(k) + "\t");
				
				sb.append(String.format("%s%n", combin.get(combin.size()-1)));
				
				bw.write(sb.toString());
			}
		} 
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		finally {
			try {
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
	}//end of printCombinations method 

}//end of ParametersGrid class
