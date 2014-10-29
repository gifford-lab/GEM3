package edu.mit.csail.cgs.utils.strings;

import java.util.ArrayList;
import java.util.Collection;

public class StringUtils {

    /**
     * same thing as perl's join()
     */
    public static String join(Collection objects, String glue) {
        StringBuilder builder = new StringBuilder();
        for (Object o : objects) {
            builder.append(o);
            builder.append(glue);
        }
        builder.setLength(builder.length() - 1);
        return builder.toString();
    }
	
	public static String padString(String s, int len) { return padString(s, ' ', len); }

	public static String padString(String s, char c, int len) { 
		if(s.length() >= len) { return s; }
		StringBuilder sb = new StringBuilder(s);
		while(sb.length() < len) { 
			sb.append(c);
		}
		return sb.toString();
	}

	/**
	 * Returns a string containing String <tt>s</tt> padded by <tt>numReps</tt> 
	 * times of the String <tt>paddingPattern</tt> from the left (<tt>directionOfPadding = -1</tt>)
	 * or the right (<tt>directionOfPadding = 1</tt>)
	 * @param s String to be padded
	 * @param paddingPattern padding pattern
	 * @param numReps number of times that the padding will be performed
	 * @param directionOfPadding direction of padding. <br>
	 * <tt> Left  = -1  </tt>   <br>
	 * <tt> Right =  1  </tt>
	 * @return
	 */
	public static String padString(String s, String paddingPattern, int numReps, int directionOfPadding)
	{
		if(numReps < 0)
			throw new IllegalArgumentException("The number of repetitions must be a non negative integer.");
		
		String expandPad = "";
		for(int i = 0; i < numReps; i++) { expandPad += paddingPattern; }
		
		switch(directionOfPadding)
		{
			// left end of String
			case -1:
			{
				s = expandPad + s; break;
			
			}
		
			// right end of String
			case 1:
			{
				s += expandPad; break;
			}
		
			default:
				throw new IllegalArgumentException("Direction of padding must be either from the left end " +
												   "of String (-1) or the right (1).");
		}
		
		return s;
	}//end of padString method
	
	/**
	 * Find all occurrences (start positions) of pattern in string (single-stranded)
	 * @param str
	 * @param pattern
	 * @return a list of all matched positions (starts, i.e. pattern_string)
	 */
	public static ArrayList<Integer> findAllOccurences (String str, String pattern){
		ArrayList<Integer> pos = new ArrayList<Integer>();
		int len = pattern.length();
		if (len > 0) {  
			int start = str.indexOf(pattern);
			while (start != -1) {
				pos.add(start);
				start = str.indexOf(pattern, start+len);
			}
		}
		return pos;
  }
}