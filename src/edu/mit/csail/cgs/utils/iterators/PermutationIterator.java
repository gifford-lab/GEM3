/*
 * Created on Mar 5, 2006
 */
package edu.mit.csail.cgs.utils.iterators;

import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class PermutationIterator implements Iterator<String> {
    
    public static void main(String[] args) { 
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        System.out.print("> "); System.out.flush();
        try {
            String line = br.readLine().trim();
            String[] array = line.split("\\s+");
            int length = Integer.parseInt(array[0]);
            int ones = Integer.parseInt(array[1]);
            System.out.println("[" + length + "," + ones + "]");
            PermutationIterator itr = new PermutationIterator(length, ones);
            while(itr.hasNext()) { 
                System.out.println(itr.next());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }
    
    private int[] array;
    private int numOnes;
    private String nextString;

    public PermutationIterator(int length, int ones) {
        array = new int[length];
        numOnes = ones;
        for(int i = 0; i < array.length; i++) { array[i] = 0; }
        for(int i = 0; i < ones; i++) { array[i] = 1; }
        nextString = getString();
        System.out.println("start: " + getString());
    }
    
    private String getString() { 
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < array.length; i++) { sb.append(String.valueOf(array[i])); }
        return sb.toString();
    }
    
    private boolean ended() { 
        for(int i = array.length-1; i > array.length-1-numOnes; i--) { 
            if(array[i] == 0) { return false; }
        }
        return true;
    }
    
    private void update() {
        int i;
        int seen = 0;
        for(i = array.length-1; i >= 0; i--) { 
            if(array[i] == 1) { 
                if(i < array.length-1 && array[i+1] == 0) { 
                    array[i+1] = 1;
                    array[i] = 0;
                    
                    for(int j = i+2; j < array.length; j++) {
                        if(j-(i+2) < seen) { 
                            array[j] = 1;
                        } else { 
                            array[j] = 0;
                        }
                    }
                    return;
                } else { 
                    seen++;
                }
            }
        }    
        
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return nextString != null;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public String next() {
        String ret = getString();
        if(ended()) { 
            nextString = null; 
        } else {
            update();
            nextString = getString();
        }
        return ret;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

}
