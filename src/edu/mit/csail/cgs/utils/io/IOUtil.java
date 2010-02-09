package edu.mit.csail.cgs.utils.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

import java.util.ArrayList;

import edu.mit.csail.cgs.utils.Utils;

/**
 * This class will hold IO operations for various types of variables: arrays, matrices, lists etc
 * @author geopapa
 *
 */
public class IOUtil {
	
	/**********************************
	 **       I OPERATIONS           **
	 **********************************/
	
  	 /////////////////////////
	 //   READ FROM FILE   //
    ////////////////////////

	public static String[] readFile2Array(String f) {
		return readFile2Array(f, "\t");
	}
	
	public static String[] readFile2Array(String f, String delim) {
		String[][] a = readFile(f, delim, -1);
		return a[0];
	}
	
	public static String[][] readFile(String f, int num_cols) {
		return readFile(f, "\t", num_cols);
	}
	
	public static String[][] readFile(String f, String delim) {
		return readFile(f, delim, 0);
	}

	/**
	 * 
	 * @param f file from where the array will be read
	 * @param delim delimiter used
	 * @param num_cols number of columns of the array
	 * @return
	 */
	public static String[][] readFile(String f, String delim, int num_cols) {
		String[][] a = new String[0][];
		try { a = readStream(new FileInputStream(new File(f)), delim, num_cols); } 
		catch (FileNotFoundException e) { e.printStackTrace(); }
		return a;
	}

	
	 ///////////////////////////
	 //   READ FROM STREAM   //
    //////////////////////////
	
	public static String[] readStream2Array(InputStream os) {
		return readStream2Array(os, "\t");
	}
	
	public static String[] readStream2Array(InputStream os, String delim) {
		String[][] a = readStream(os, delim, -1);
		return a[0];
	}
	
	public static String[][] readStream(InputStream os, int num_cols) {
		return readStream(os, "\t", num_cols);
	}
	
	public static String[][] readStream(InputStream os, String delim) {
		return readStream(os, delim, 0);
	}

	/**
	 * 
	 * @param os The output stream where the array will be stored
	 * @param a array to be stored
	 * @param delim delimiter used
	 * @param num_cols number of array elements in each line (row) of the stream
	 */
	public static String[][] readStream(InputStream os, String delim, int num_cols) {
		String[][] a = new String[0][];
		ArrayList<String[]> al = new ArrayList<String[]>();
		BufferedReader br = null;
		
		try {
			
			InputStreamReader osr = new InputStreamReader(os);
			br = new BufferedReader(osr);
			
			if( num_cols > 0 ) {		
				String str;
				int count = 0;
				String[] leftover = new String[0];
				while((str = br.readLine()) != null) {
					String[] tokens = str.split(delim);
					String[] els = concat_strarrays(leftover, tokens); 
					count = els.length;
					
					int quot = count/num_cols;
					
					int ind = 0;
					if(quot != 0) {
						for(int i = 0; i < quot; i++) { 
							String[] curr_els = new String[num_cols];
							for(int j = 0; j < num_cols; j++) { curr_els[j] = els[ind++]; }
							al.add(curr_els);
						}
						count %= num_cols;
					}
					
					leftover = new String[count];
					for(int k = 0; k < count; k++)
						leftover[k] = els[ind++];	
				}
				al.add(leftover);
			}
			
			else if( num_cols == 0 ) {
				String str;
				while((str = br.readLine()) != null) {
					String[] tokens = str.split(delim);
					al.add(tokens);
				}
			}
			
			else if( num_cols == -1 ) {
				ArrayList<String> temp_list = new ArrayList<String>();
				String str;
				while((str = br.readLine()) != null) {
					String[] tokens = str.split(delim);
					for(int i = 0; i < tokens.length; i++) { temp_list.add(tokens[i]); }
				}
				al.add(temp_list.toArray(new String[0]));
			}
			else { throw new IllegalArgumentException("Invalid value for num_cols. The valid values are -1 (one row), 0 (as they are on the file), k (columns in each row)."); }
	
			a = new String[al.size()][];
			for(int i = 0; i < a.length; i++) { a[i] = al.get(i); }
		}
		catch (IOException e) { e.printStackTrace(); }
		finally {
			try { br.close(); } 
			catch (IOException e) { e.printStackTrace(); }
		}
		
		return a;
	}//end of readStream method
	
	
	private static String[] concat_strarrays(String[] ar1, String[] ar2) {
		String[] ar = new String[ar1.length+ar2.length];
		int count = 0;
		for(int i = 0; i < ar1.length; i++) { ar[count++] = ar1[i]; }
		for(int i = 0; i < ar2.length; i++) { ar[count++] = ar2[i]; }
		return ar;
	}

	
	
	/**********************************
	 **       O OPERATIONS           **
	 **********************************/
	
   	 ///////////////////////
	 //   WRITE 2 FILE   //
    ///////////////////////
	public static <T> void write2file(String f, T[][][] a) {
		write2file(f, a, "\t");
	}//end of write2file method
	
	public static <T> void write2file(String f, T[][][] a, String delim) {
		try { write2stream(new FileOutputStream(new File(f)), a, delim); } 
		catch (FileNotFoundException e) { e.printStackTrace(); }				
	}//end of write2file method
	
	public static <T> void write2file(String f, T[][] a) {
		write2file(f, a, "\t");
	}//end of write2file method
	
	public static <T> void write2file(String f, T[][] a, String delim) {
		try { write2stream(new FileOutputStream(new File(f)), a, delim); } 
		catch (FileNotFoundException e) { e.printStackTrace(); }		
	}//end of write2file method
	
	public static <T> void write2file(String f, T[] a) {
		write2file(f, a, "\t", -1);	
	}//end of write2file method	
	
	public static <T> void write2file(String f, T[] a, int num_cols) {
		write2file(f, a, "\t", num_cols);
	}//end of write2file method
	
	public static <T> void write2file(String f, T[] a, String delim) {
		write2file(f, a, delim, -1);
	}//end of write2file method
	
	public static <T> void write2file(String f, T[] a, String delim, int num_cols) {
		try { write2stream(new FileOutputStream(new File(f)), a, delim, num_cols); } 
		catch (FileNotFoundException e) { e.printStackTrace(); }
	}//end of write2file method
	
	
  	  ///////////////////////
	 //   WRITE 2 STREAM   //
    ///////////////////////
	public static <T> void write2stream(OutputStream os, T[][][] a) {
		write2stream(os, a, "\t");
	}//end of write2stream method
	
	public static <T> void write2stream(OutputStream os, T[][][] a, String delim) {
		BufferedWriter bw = null;
		try {
			OutputStreamWriter osw = new OutputStreamWriter(os);
			bw = new BufferedWriter(osw);
			
			for(int k = 0; k < a.length; k++) {
				for(int i = 0; i < a[k].length; i++) {
					for(int j = 0; j < a[k][i].length-1; j++)
						bw.write(a[k][i][j].toString() + delim);
					
					bw.write(a[k][i][a[k][i].length-1].toString() + "\n");
				}
				if( k < a.length-1) { bw.write("\n"); }
			}
			if(a.length == 0) { bw.write("\n"); }
		}
		catch (IOException e) { e.printStackTrace(); }
		finally {
			try {
				bw.close();
			} 
			catch (IOException e) { e.printStackTrace(); }
		}
	}//end of write2stream method
	
	public static <T> void write2stream(OutputStream os, T[][] a) {
		write2stream(os, a, "\t");
	}//end of write2stream method	
	
	public static <T> void write2stream(OutputStream os, T[][] a, String delim) {
		BufferedWriter bw = null;
		try {
			OutputStreamWriter osw = new OutputStreamWriter(os);
			bw = new BufferedWriter(osw);
			
			for(int i = 0; i < a.length; i++) {
				for(int j = 0; j < a[i].length-1; j++)
					bw.write(a[i][j].toString() + delim);
				
				bw.write(a[i][a[i].length-1].toString() + "\n");
			}
			if(a.length == 0) { bw.write("\n"); }
		}
		catch (IOException e) { e.printStackTrace(); }
		finally {
			try {
				bw.close();
			} 
			catch (IOException e) { e.printStackTrace(); }
		}
	}//end of write2stream method
	
	public static <T> void write2stream(OutputStream os, T[] a) {
		write2stream(os, a, "\t", -1);	
	}//end of write2stream method	
	
	public static <T> void write2stream(OutputStream os, T[] a, int num_cols) {
		write2stream(os, a, "\t", num_cols);
	}//end of write2stream method
	
	public static <T> void write2stream(OutputStream os, T[] a, String delim) {
		write2stream(os, a, delim, -1);
	}//end of write2stream method
	
	/**
	 * 
	 * @param <T>
	 * @param os The output stream where the array will be stored
	 * @param a array to be stored
	 * @param delim demimiter used
	 * @param num_cols number of array elements in each line (row) of the stream
	 */
	public static <T> void write2stream(OutputStream os, T[] a, String delim, int num_cols) {
		BufferedWriter bw = null;
		try {
			OutputStreamWriter osw = new OutputStreamWriter(os);
			bw = new BufferedWriter(osw);
			
			if( num_cols > 0) {
				for(int n = 0; n < a.length-1; n++) {
					if((n+1)%num_cols == 0) { bw.write(a[n].toString() + "\n"); }
					else 				    { bw.write(a[n].toString() + delim); }
				}				
			}
			else {
				for(int n = 0; n < a.length-1; n++) {
					bw.write(a[n].toString() + delim);
				}								
			}
			
			if(a.length > 0) { bw.write(a[a.length-1].toString() + "\n"); }
			else             { bw.write("\n");                            }
		}
		catch (IOException e) { e.printStackTrace(); }
		finally {
			try {
				bw.close();
			} 
			catch (IOException e) { e.printStackTrace(); }
		}
	}//end of write2stream method


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		String f = "/Users/gio_fou/Desktop/skatotest.txt";
		
		String[][] rr = readFile(f, ",", -1);
		
		int[][][] a = {
		             { {3, 2, 1, 5, 3, 6, 2, 7}, 
				       {1, 6, 3, 7, 2, 1} },
			           { {-2, -4, 2, 1, 5}, 
					     {0, 4, 2, -1} },				       
		            };
		String delim = ",";
		
		int[][] b = { {3, 1, 2, 5}, {4, -1, 4} };
		
		write2file(f, Utils.prim2ref(a), delim);
		write2stream(System.out, Utils.prim2ref(b));
		
	}

}
