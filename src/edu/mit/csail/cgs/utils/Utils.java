package edu.mit.csail.cgs.utils;

import java.util.*;
import java.io.*;

public abstract class Utils {
    
    public static Map<String,String> readMap(String filename) 
    throws FileNotFoundException, IOException {
        return readMap(filename,false);
    }
    
    public static Map<String,String> readMap(String filename, boolean reverse) 
    throws FileNotFoundException, IOException {
        BufferedReader r = new BufferedReader(new FileReader(filename));
        String line;
        HashMap<String,String> map = new HashMap<String,String>();
        
        while ((line = r.readLine()) != null) {
            StringTokenizer tok = new StringTokenizer(line,"\t");
            String key,value;
            if (tok.countTokens() >= 2) {
                key = tok.nextToken();
                value = tok.nextToken();
                if (reverse) {
                    map.put(value,key);
                } else {
                    map.put(key,value);
                }
            }
        }
        r.close();
        return map;
    }
    
    /**
     * Saves a map to a file with the given name -- by default, calls the corresponding 
     * saveMap() function with a "reverse" parameter of <i>false</i>.
     * 
     * @param filename The name of the file to which to save.
     * @param map The String-to-String map to save.
     * @throws IOException
     */
    public static void saveMap(String filename, Map<String,String> map) 
    throws IOException {
        saveMap(filename,map,false);
    }
    
    /**
     * Saves a map to a file with the given name.
     * 
     * @param filename The name of the file to which to save.
     * @param map The String-to-String map to save.
     * @param reverse If true, saves a map not of key:value, but of value:key.
     * @throws IOException
     */
    public static void saveMap(String filename, Map<String,String> map, 
            boolean reverse) 
    throws IOException {
        
        BufferedWriter w = new BufferedWriter(new FileWriter(filename));
        for(String k : map.keySet()) { 
            String v = map.get(k);
            if (reverse) {
                w.write(v + "\t" + k + "\n");
            } else {
                w.write(k + "\t" + v + "\n");
            }
        }
        w.close();
    }
    
    public static String transformDuration(long msec_duration) { return transformDuration(msec_duration, "dhms");   }
    
    public static String transformDuration(long msec_duration, String format) {
     	if( msec_duration < 0) { throw new IllegalArgumentException("msec_duration must be a non-negative integer."); }
     	
     	if(!format.matches("^[DdHhMmSs]?[DdHhMmSs]?[DdHhMmSs]?[DdHhMmSs]?$")) 
     		throw new IllegalArgumentException("The format has to be a set of: dhms. " +
     										   "d: for days, h: for hours, m: for minutes, s: for seconds." +
     				                           "E.g.: dh, ds, hms and so on.");
     	
     	StringBuilder sb = new StringBuilder();
     	long quot;
     	
     	long sec_duration_long   = msec_duration/1000;
     	double sec_duration      = msec_duration/1000.0;
     	sec_duration            -= sec_duration_long;
     	
     	// Get day information
     	if(format.matches("[Dd]+")) { return String.format("%.1f d", (double)sec_duration_long/86400); }
     	else if(format.matches(".*[Dd]+.*")) {
     		if( (quot = sec_duration_long/86400) != 0) { 
     			sb.append(quot + " d, "); 
     			sec_duration_long -= quot*86400;
     		}
     		format = format.replaceAll("[Dd]", "");
     	}
     	
     	// Get hour information
     	if(format.matches("[Hh]+")) { return sb.append(String.format("%.1f h", (double)sec_duration_long/3600)).toString(); }
     	else if(format.matches(".*[Hh]+.*")) {
     		if( (quot = sec_duration_long/3600) != 0) { 
     			sb.append(quot + " h, "); 
     			sec_duration_long -= quot*3600;
     		}
     		format = format.replaceAll("[Hh]", "");
     	}
     	
     	// Get minute information
     	if(format.matches("[Mm]+")) { return sb.append(String.format("%.1f m", (double)sec_duration_long/60)).toString(); }
     	else if(format.matches(".*[Mm]+.*")) {
     		if( (quot = sec_duration_long/60) != 0) {
     			sb.append(quot + " m, ");
     			sec_duration_long -= quot*60;
     		}
     		format = format.replaceAll("[Mm]", "");
     	}
     	
     	sec_duration += sec_duration_long;
     	
     	// Get second information
     	sb.append(String.format("%.2f s", sec_duration));

     	return sb.toString();
     }//end of transformDuration method
    
    
    /******************************************
     **      PRIMITIVE 2 REFERENCE TYPE      **
     ******************************************/
    
    // Byte
    public static Byte[][][] prim2ref(byte[][][] a) {
    	Byte[][][] b = new Byte[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Byte[][] prim2ref(byte[][] a) {
    	Byte[][] b = new Byte[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method
    
    public static Byte[] prim2ref(byte[] a) {
    	Byte[] b = new Byte[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Short
    public static Short[][][] prim2ref(short[][][] a) {
    	Short[][][] b = new Short[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Short[][] prim2ref(short[][] a) {
    	Short[][] b = new Short[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method
    
    public static Short[] prim2ref(short[] a) {
    	Short[] b = new Short[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Integer
    public static Integer[][][] prim2ref(int[][][] a) {
    	Integer[][][] b = new Integer[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Integer[][] prim2ref(int[][] a) {
    	Integer[][] b = new Integer[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method
    
    public static Integer[] prim2ref(int[] a) {
    	Integer[] b = new Integer[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Long
    public static Long[][][] prim2ref(long[][][] a) {
    	Long[][][] b = new Long[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Long[][] prim2ref(long[][] a) {
    	Long[][] b = new Long[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method
    
    public static Long[] prim2ref(long[] a) {
    	Long[] b = new Long[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Float
    public static Float[][][] prim2ref(float[][][] a) {
    	Float[][][] b = new Float[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Float[][] prim2ref(float[][] a) {
    	Float[][] b = new Float[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method

    public static Float[] prim2ref(float[] a) {
    	Float[] b = new Float[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Double
    public static Double[][][] prim2ref(double[][][] a) {
    	Double[][][] b = new Double[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Double[][] prim2ref(double[][] a) {
    	Double[][] b = new Double[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method    

    public static Double[] prim2ref(double[] a) {
    	Double[] b = new Double[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    // Character
    public static Character[][][] prim2ref(char[][][] a) {
    	Character[][][] b = new Character[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = prim2ref(a[k]);
    	return b;
    }//end of prim2ref method

    public static Character[][] prim2ref(char[][] a) {
    	Character[][] b = new Character[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = prim2ref(a[t]);
    	return b;
    }//end of prim2ref method        

    public static Character[] prim2ref(char[] a) {
    	Character[] b = new Character[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of prim2ref method
    
    
    
    /******************************************
     **      REFERENCE 2 PRIMITIVE TYPE      **
     ******************************************/
    
    // byte
    public static byte[][][] ref2prim(Byte[][][] a) {
    	byte[][][] b = new byte[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static byte[][] ref2prim(Byte[][] a) {
    	byte[][] b = new byte[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method
    
    public static byte[] ref2prim(Byte[] a) {
    	byte[] b = new byte[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // short
    public static short[][][] ref2prim(Short[][][] a) {
    	short[][][] b = new short[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static short[][] ref2prim(Short[][] a) {
    	short[][] b = new short[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method
    
    public static short[] ref2prim(Short[] a) {
    	short[] b = new short[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // int
    public static int[][][] ref2prim(Integer[][][] a) {
    	int[][][] b = new int[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static int[][] ref2prim(Integer[][] a) {
    	int[][] b = new int[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method
    
    public static int[] ref2prim(Integer[] a) {
    	int[] b = new int[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // long
    public static long[][][] ref2prim(Long[][][] a) {
    	long[][][] b = new long[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static long[][] ref2prim(Long[][] a) {
    	long[][] b = new long[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method
    
    public static long[] ref2prim(Long[] a) {
    	long[] b = new long[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // float
    public static float[][][] ref2prim(Float[][][] a) {
    	float[][][] b = new float[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static float[][] ref2prim(Float[][] a) {
    	float[][] b = new float[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method

    public static float[] ref2prim(Float[] a) {
    	float[] b = new float[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // double
    public static double[][][] ref2prim(Double[][][] a) {
    	double[][][] b = new double[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static double[][] ref2prim(Double[][] a) {
    	double[][] b = new double[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method    

    public static double[] ref2prim(Double[] a) {
    	double[] b = new double[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
    
    
    // char
    public static char[][][] ref2prim(Character[][][] a) {
    	char[][][] b = new char[a.length][][];
    	for(int k = 0; k < b.length; k++)
    		b[k] = ref2prim(a[k]);
    	return b;
    }//end of ref2prim method

    public static char[][] ref2prim(Character[][] a) {
    	char[][] b = new char[a.length][];
    	for(int t = 0; t < b.length; t++)
    		b[t] = ref2prim(a[t]);
    	return b;
    }//end of ref2prim method        

    public static char[] ref2prim(Character[] a) {
    	char[] b = new char[a.length];
    	for(int i = 0; i < b.length; i++) { b[i] = a[i]; }
    	return b;
    }//end of ref2prim method
   
}// end of Utils class
