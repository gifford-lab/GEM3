/*
 * Created on Oct 14, 2005
 */
package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.text.*;

/**
 * @author tdanford
 */
public class PairwiseLetterDistribution {
    
    private static NumberFormat nf;
    private static int numberLength;
    
    static { 
        nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(3);
        nf.setMinimumFractionDigits(3);
        nf.setMinimumIntegerDigits(1);
        numberLength = 5;
    }

    private int totalCount;
    private int[][] countArray;
    private int[][] marginalArray;
    private char[] letters;
    private Map<Character,Integer> indexMap;

    public PairwiseLetterDistribution(char[] c) {
        letters = new char[c.length];
        indexMap = new TreeMap<Character,Integer>();
        countArray = new int[c.length][c.length];
        marginalArray = new int[2][c.length];
        for(int i = 0; i < c.length; i++) {
            letters[i] = Character.toUpperCase(c[i]);
            indexMap.put(Character.toUpperCase(c[i]), i);
            marginalArray[0][i] = 0; 
            marginalArray[1][i] = 0;
            for(int j = 0; j < c.length; j++) { 
                countArray[i][j] = 0;
            }
        }
        totalCount = 0;
    }
    
    public int getConservedCount() {
        int sum = 0;
        for(int i = 0; i < letters.length; i++) { 
            sum += countArray[i][i];
        }
        
        return sum;
    }
    
    public double getConservedFraction() { 
        if(totalCount == 0) { return 0.0; }
        int conservedCount = getConservedCount();
        return (double)conservedCount / (double)totalCount;
    }
    
    public double getInformation() { 
        double sum = 0.0;
        for(int i = 0; i < letters.length; i++) { 
            for(int j = 0; j < letters.length; j++) { 
                double p = getFraction(letters[i], letters[j]);
                if(p != 0.0) { 
                    sum += (p * Math.log(p));
                }
            }
        }
        return -sum;
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        
        sb.append("# Count: " + totalCount + "\n");
        sb.append("Info (Bits): " + getInformation() + "\n");
        sb.append("Conservation %: " + getConservedFraction());
        
        sb.append(spaceBlock(2)); sb.append("|");
        int hw = numberLength/2;
        for(int i = 0; i < letters.length; i++) { 
            sb.append(spaceBlock(hw));
            sb.append(letters[i]);
            sb.append(spaceBlock(numberLength-hw));            
        }
        sb.append("\n");
        
        sb.append("---");
        for(int i = 0; i < letters.length; i++) { 
            sb.append(block(numberLength+1, '-'));
        }
        sb.append("\n");
        
        for(int i = 0; i < letters.length; i++) { 
            sb.append(letters[i] + " |");
            for(int j = 0; j < letters.length; j++) { 
                sb.append(nf.format(getFraction(letters[i], letters[j])));
                sb.append(" ");
            }
            sb.append(" : ");
            sb.append(nf.format(getMarginalFraction(letters[i], true)));
            sb.append("\n");
        }
        
        sb.append(spaceBlock(3));
        for(int j = 0; j < letters.length; j++) { 
            sb.append(nf.format(getMarginalFraction(letters[j], false)));
            sb.append(" ");
        }
        sb.append("\n");
        
        return sb.toString();
    }
    
    private String spaceBlock(int len) { return block(len, ' '); }
    
    private String block(int len, char c) { 
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < len; i++) { sb.append(c); }
        return sb.toString();
    }
    
    public int getTotalCount() { return totalCount; }
    public boolean hasChar(char c) { return indexMap.containsKey(Character.toUpperCase(c)); }
    public int getCharIndex(char c) { return indexMap.get(Character.toUpperCase(c)); }
    
    public int getCount(char c1, char c2) { 
        int i1 = getCharIndex(c1), i2 = getCharIndex(c2);
        return countArray[i1][i2];
    }
    
    public void addCount(char c1, char c2, int c) { 
		int i1 = getCharIndex(c1), i2 = getCharIndex(c2);
		//System.out.print(c1 + "," + c2 + " " + countArray[i1][i2] + " --> ");
        countArray[i1][i2] += c;
		marginalArray[0][i1] += c;
		marginalArray[1][i2] += c;
		totalCount += c;
		//System.out.println(countArray[i1][i2]);
    }
    
    public double getFraction(char c1, char c2) { 
        if(totalCount == 0) { return 0.0; }
        return (double)getCount(c1, c2) / (double)totalCount;
    }
    
    public int getMarginalCount(char c, boolean first) { 
        return marginalArray[first ? 0 : 1][getCharIndex(c)];
    }
    
    public double getMarginalFraction(char c, boolean first) { 
        if(totalCount == 0) { return 0.0; }
        return (double)getMarginalCount(c, first) / (double)totalCount;
    }
}
