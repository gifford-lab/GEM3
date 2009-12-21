/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import java.util.*;
import java.util.regex.*;
import java.lang.reflect.*;
import java.io.*;

import edu.mit.csail.cgs.utils.Accumulator;
import edu.mit.csail.cgs.utils.Function;
import edu.mit.csail.cgs.utils.PackedBitVector;
import edu.mit.csail.cgs.utils.Predicate;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.ModelFieldAnalysis;
import edu.mit.csail.cgs.utils.models.ModelInput;
import edu.mit.csail.cgs.utils.models.ModelInputIterator;
import edu.mit.csail.cgs.utils.models.ModelOutput;
import edu.mit.csail.cgs.utils.models.ModelInput.LineReader;
import edu.mit.csail.cgs.utils.models.ModelOutput.LineWriter;

/**
 * DataFrame holds a collection of Model objects -- intuitively, a DataFrame is a matrix 
 * of data, although the matrix data is not necessarily numerical.  Each internal Model is a
 * row in the matrix, and each field in the Model is a column.  The DataFrame gives accessor
 * methods for getting data by either row or column, for reading and writing tables to files,
 * for assembling matrices of numerical data, etc. 
 * 
 * @author tdanford
 * @param <T>  The particular subclass of Model whose objects make up the rows of the frame.
 */
public class DataFrame<T extends Model> {

	private Class<T> cls;
	private ArrayList<T> objects;
	private File file;
	private ModelFieldAnalysis<T> fieldAnalysis;
	
	/**
	 * Reads the frame data from a particular file.  
	 * 
	 * @param cls
	 * @param file
	 * @throws IOException
	 */
	public DataFrame(Class<T> cls, File file) throws IOException { 
		this.cls = cls;
		this.file = file;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
		objects = parse(file, true);
	}
	
	/**
	 * Reads the frame data from a particular file.  
	 * 
	 * @param cls
	 * @param file
	 * @throws IOException
	 */
	public DataFrame(Class<T> cls, File file, String... cols) throws IOException { 
		this.cls = cls;
		this.file = file;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
		objects = parse(file, true, cols);
	}

	public DataFrame(Class<T> cls, File file, boolean header, String... cols) throws IOException { 
		this.cls = cls;
		this.file = file;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
		objects = parse(file, header, cols);
	}
	
	/**
	 * Creates a DataFrame object from the given iterator over values.
	 * 
	 * @param cls
	 * @param objs
	 */
	public DataFrame(Class<T> cls, Iterator<T> objs) { 
		this.cls = cls;
		objects = new ArrayList<T>();
		while(objs.hasNext()) { 
			objects.add(objs.next());
		}
		this.file = null;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
	}
	
	/**
	 * Creates a DataFrame object from the given collection of values.
	 * 
	 * @param cls
	 * @param objs
	 */
	public DataFrame(Class<T> cls, Collection<T> objs) { 
		this.cls = cls;
		objects = new ArrayList<T>(objs);
		this.file = null;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
	}
	
	/**
	 * Creates an empty DataFrame object.  
	 * 
	 * @param cls
	 */
	public DataFrame(Class<T> cls) { 
		this.cls = cls;
		objects = new ArrayList<T>();
		file = null;
		fieldAnalysis = new ModelFieldAnalysis<T>(cls);
	}
	
	public void loadJSON(InputStream is) { 
		ModelInput<T> input = new ModelInput.LineReader<T>(cls, is);
		ModelInputIterator<T> itr = new ModelInputIterator<T>(input);
		addObjects(itr);
	}
	
	public void saveJSON(OutputStream os) { 
		ModelOutput<T> output = new ModelOutput.LineWriter<T>(os);
		for(T val : objects) { output.writeModel(val); }
	}
	
	public Iterator<T> iterator() { return objects.iterator(); }
	
	public PackedBitVector getMask(Predicate<T> pred) { 
		PackedBitVector pbv = new PackedBitVector(objects.size());
		for(int i = 0; i < objects.size(); i++) { 
			if(pred.accepts(objects.get(i))) { 
				pbv.turnOnBit(i);
			}
		}
		return pbv;
	}
	
	/**
	 * Removes any item that satisfies the given predicate from this DataFrame.  
	 * Extracted items are collated in a separate DataFrame<T>, which is then
	 * returned from this method.
	 * 
	 * @param pred The indicator for which elements are to be removed.
	 * @return A DataFrame of the removed elements.  
	 */
	public DataFrame<T> extract(Predicate<T> pred) { 
		DataFrame<T> extracted = new DataFrame<T>(cls);
		Iterator<T> itr = objects.iterator();
		while(itr.hasNext()) { 
			T value = itr.next();
			if(pred.accepts(value)) { 
				itr.remove();
				extracted.addObject(value);
			}
		}
		return extracted;
	}
	
	/**
	 * A Transformation<T,S> object turns objects of type T into objects of type S.  
	 * The transform() method takes a transformation, where T is the model-type of 
	 * this DataFrame, and returns a new DataFrame of S objects, where each object
	 * corresponds to a transformed version of the original T from this frame.  
	 * 
	 * @param <S>
	 * @param trans
	 * @return
	 */
	public <S extends Model> DataFrame<S> transform(Transformation<T,S> trans) { 
		DataFrame<S> df = new DataFrame<S>(trans.toClass());
		for(T val : objects) { 
			df.addObject(trans.transform(val));
		}
		return df;
	}
	
	public DataFrame<T> extend(DataFrame<T> df) { 
		objects.addAll(df.objects);
		return this;
	}
	
	public <R extends Model, S extends Model> 
		DataFrame<R> join(Class<R> rClass, 
				DataFrame<S> outerFrame, 
				String fieldName, 
				String innerName, 
				String outerName) { 
		
		ModelFieldAnalysis<S> analysisS = outerFrame.fieldAnalysis;
		Field fieldS = analysisS.findField(fieldName);

		ModelFieldAnalysis<T> analysisT = fieldAnalysis;
		Field fieldT = analysisT.findField(fieldName);
		
		ModelFieldAnalysis<R> analysisR = new ModelFieldAnalysis<R>(rClass);
		Field rJoin = analysisR.findField(fieldName);
		Field rInner = analysisR.findField(innerName);
		Field rOuter = analysisR.findField(outerName);
		
		if(fieldS == null || fieldT == null || rJoin == null) { 
			throw new IllegalArgumentException(fieldName);
		}
		
		LinkedList<R> joinedValues = new LinkedList<R>();
		HashMap<Object,ArrayList<S>> hashS = new HashMap<Object,ArrayList<S>>(); 
		for(int i = 0; i < outerFrame.size(); i++) { 
			S outerValue = outerFrame.object(i);
			try {
				Object key = fieldS.get(outerValue);
				
				if(!hashS.containsKey(key)) { 
					hashS.put(key, new ArrayList<S>());
				}
				hashS.get(key).add(outerValue);
				
			} catch (IllegalAccessException e) {
				e.printStackTrace();
				throw new IllegalArgumentException(String.format("Field %s " +
						"was illegally accessed in Model %s", fieldName,
						outerFrame.getModelClass().getSimpleName()));
			}
		}
		
		for(int i = 0; i < size(); i++) { 
			T innerValue = object(i);
			try {
				Object key = fieldT.get(innerValue);
				if(hashS.containsKey(key)) { 
					for(S outerValue : hashS.get(key)) { 
						R joined = (R)rClass.newInstance();
						
						rJoin.set(joined, key);
						
						if(rOuter != null) { 
							rOuter.set(joined, outerValue);
						}
						
						if(rInner != null) { 
							rInner.set(joined, innerValue);
						}
						
						joinedValues.add(joined);
					}
				}
			} catch (IllegalAccessException e) {
				e.printStackTrace();
				throw new IllegalArgumentException(String.format("Field %s " +
						"was illegally accessed in Model %s", fieldName,
						getModelClass().getSimpleName()));
			} catch (InstantiationException e) {
				e.printStackTrace();
				throw new IllegalArgumentException(String.format(
						"Couldn't instantiate Model class %s", 
						rClass.getSimpleName()));
			}
		}

		DataFrame<R> joinFrame = new DataFrame<R>(rClass, joinedValues);
		return joinFrame; 
	}
	
	public DataFrame<T> filter(Predicate<T> p) { 
		DataFrame<T> df = new DataFrame<T>(cls);
		for(T val : objects) {  
			if(p.accepts(val)) { 
				df.addObject(val);
			}
		}
		return df;
	}
	
	public void apply(Accumulator<T> acc) { 
		for(T val : objects) { 
			acc.accumulate(val);
		}
	}
	
	/**
	 * Saves the data back out to a particular file -- this file then becomes the DataFrame's
	 * "default" file, and is the target when future calls are made to the save() method 
	 * (no parameters).  
	 * 
	 * @param f
	 * @throws IOException
	 */
	public void save(File f) throws IOException {
		writeTable(objects, fieldAnalysis.getFieldNames(), f);
		file = f;
	}
	
	/**
	 * Saves the data back out to the "default" file (either the file from which the data
	 * was loaded, or the file that was given as a parameter to the last call to save(File).)
	 * 
	 * @throws IOException
	 */
	public void save() throws IOException { 
		save(file);
	}

	/**
	 * Adds new rows to the DataFrame object.
	 *  
	 * @param objs
	 */
	public void addObjects(Iterator<T> objs) { 
		while(objs.hasNext()) { 
			addObject(objs.next());
		}
	}
	
	public File getFile() { return file; }
	public Class<T> getModelClass() { return cls; }
	
	public void addObjects(Collection<T> objs) { objects.addAll(objs); }
	public void addObject(T obj) { objects.add(obj); }
	
	public Vector<String> getFields() { return fieldAnalysis.getFieldNames(); }
	public T object(int i) { return objects.get(i); }
	public int size() { return objects.size(); }
	
	/**
	 * For a particular named field in the model class, this returns a Set of all distinct
	 * values of that field among the rows of this DataFrame.
	 *   
	 * @param <FT>
	 * @param fieldName
	 * @return
	 */
	public <FT> Set<FT> fieldValues(String fieldName) {  
		HashSet<FT> values = new HashSet<FT>();
		try {
			Field field = cls.getField(fieldName);
			for(T obj : objects) {
				FT objFieldValue = (FT)field.get(obj);
				values.add(objFieldValue);
			}
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(fieldName);
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return values; 
	}
	
	public Double pearsonCorrelation(String xCol, String yCol) { 
		Double[] x = asVector(xCol);
		Double[] y = asVector(yCol);
		
		double xSum = 0.0, ySum = 0.0, xSqSum = 0.0, ySqSum = 0.0, xySum = 0.0;
		int N = 0;
		
		for(int i = 0; i < x.length; i++) { 
			if(x[i] != null && y[i] != null) { 
				N += 1;
				xSum += x[i];
				ySum += y[i];
				xSqSum += (x[i] * x[i]);
				ySqSum += (y[i] * y[i]);
				xySum += (x[i] * y[i]);
			}
		}
		
		/**
		 *                        (n \sum xy - \sum x \sum y) 
		 * r_xy = -------------------------------------------------------------
		 *        sqrt(n \sum x^2 - (\sum x)^2) * sqrt(n \sum y^2 - (\sum y)^2) 
		 */
		
		if(N == 0) { throw new IllegalStateException("No values for correlation."); }
		
		double n = (double)N;
		double numer = n * xySum - (xSum * ySum);
		double denom1 = Math.sqrt(n * xSqSum - (xSum*xSum));
		double denom2 = Math.sqrt(n * ySqSum - (ySum * ySum));
		double denom = denom1 * denom2;
		
		return numer / denom;
	}
	
	public Double mean(String col) { 
		int count = 0;
		Double sum = null;
		Field field = fieldAnalysis.findField(col);
		
		if(field != null && Model.isSubclass(field.getType(), Number.class)) { 
			for(T obj : objects) { 
				try {
					Number num = (Number)field.get(obj);
					if(num != null) { 
						count += 1;
						sum = (sum == null ? num.doubleValue() : sum + num.doubleValue());
					}
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
			}
		}
		return count > 0 ? sum / (double)count : sum;
	}
	
	public Double variance(String col) { 
		Double mean = mean(col);
		if(mean == null) { return null; }
		return squaredError(col, new Function.Constant<T,Double>(mean));
	}
	
	public Double squaredError(String col, Function<T,Double> baseliner) { 
		int count = 0;
		Double sum = null;
		Field field = fieldAnalysis.findField(col);
		
		if(field != null && Model.isSubclass(field.getType(), Number.class)) { 
			for(T obj : objects) { 
				try {
					Number num = (Number)field.get(obj);
					if(num != null) { 
						count += 1;
						double val = num.doubleValue();
						double baseline = baseliner.valueAt(obj);
						double diff = val-baseline;
						double diff2 = diff*diff;
						sum = (sum == null ? diff2 : sum + diff2);
					}
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
			}
		}
		return count > 0 ? sum / (double)count : sum;
	}
	
	public Double[][] asMatrix(String... cols) { 
		return asMatrix(objects, cols);
	}
	
	public Double[] asVector(String fieldName) { 
		return asVector(objects, fieldName);
	}
	
	public Double[] asVector(ArrayList<T> rows, String fieldName) { 
		Double[] array = new Double[rows.size()];
		for(int i = 0; i < rows.size(); i++) { 
			T row = rows.get(i);
			try {
				Field field = cls.getField(fieldName);
				Class type = field.getType();
				if(Model.isSubclass(type, Number.class)) { 
					Object value = field.get(row);

					if(value == null){ 
						array[i] = null; 
					} else { 
						array[i] = ((Number)value).doubleValue();
					}
				} else { 
					array[i] = null;
				}

			} catch (NoSuchFieldException e) {
				e.printStackTrace();
				array[i] = null;
			} catch (IllegalAccessException e) {
				e.printStackTrace();
				array[i] = null;
			}
		}
		return array;
	}
	
	public Double[][] asMatrix(ArrayList<T> rows, String... cols) { 
		Vector<String> colFields = new Vector<String>();
		for(int i = 0; i < cols.length; i++) { colFields.add(cols[i]); }
		return asMatrix(rows, colFields);
	}
	
	public Double[][] asMatrix(ArrayList<T> rows, Vector<String> colFields) { 
		Double[][] array = new Double[rows.size()][colFields.size()];
		
		for(int i = 0; i < rows.size(); i++) { 
			T row = rows.get(i);
			for(int j = 0; j < colFields.size(); j++) { 
				String fieldName = colFields.get(j);
				try {
					Field field = cls.getField(fieldName);
					Class type = field.getType();
					if(Model.isSubclass(type, Number.class)) { 
						Object value = field.get(row);
						
						if(value == null){ 
							array[i][j] = null; 
						} else { 
							array[i][j] = ((Number)value).doubleValue();
						}
					} else { 
						array[i][j] = null;
					}
					
				} catch (NoSuchFieldException e) {
					e.printStackTrace();
					array[i][j] = null;
				} catch (IllegalAccessException e) {
					e.printStackTrace();
					array[i][j] = null;
				}
			}
		}
		
		return array;
	}
	
	/** Parsing Code **/
	
	private void writeTable(Collection<T> lines, Vector<String> fields, File f) throws IOException { 
		PrintStream ps = new PrintStream(new FileOutputStream(f));
		for(int i = 0; i < fields.size(); i++) { 
			if(i > 0) { ps.print("\t"); } 
			ps.print(fields.get(i));
		}
		ps.println();
		
		for(T line : lines) { 
			writeLine(line, fields, ps);
		}
		ps.close();
	}
	
	private ArrayList<T> parse(File f, boolean header) throws IOException { 
		return parse(f, header, (String[])null);
	}
	
	private ArrayList<T> parse(File f, boolean header, String... fieldArray) throws IOException {
		ArrayList<T> values = new ArrayList<T>();
		BufferedReader br = new BufferedReader(new FileReader(f));

		Vector<String> fields = new Vector<String>();

		String line = null;
		String[] array = null;
		String sep = "\\s+";
		
		if(header) { 
			line = br.readLine();
		}
		
		if(fieldArray != null && fieldArray.length > 0) {
			for(int i = 0; i < fieldArray.length; i++) { 
				fields.add(fieldArray[i]);
			}
		} else if (header && line != null) {  
			array = line.split(sep);
			for(int i = 0; i < array.length; i++) { fields.add(array[i]); }
		}

		Vector<Boolean> quoted = new Vector<Boolean>();
		for(String fo : fields) { 
			boolean isquoted = fieldAnalysis.getStaticSwitch(String.format("quote_%s", fo), false);
			quoted.add(!isquoted);
		}
		
		int ignored = 0;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			
			if(line.length() > 0) {
				array = line.split(sep);
				T value = parseLine(array, fields, quoted);
				
				if(value != null) { 
					values.add(value); 
				} else { 
					ignored += 1;
				}
			}
		}
		
		System.out.println(String.format("Parsed %d lines from %s", values.size(), f.getName()));
		if(ignored > 0) { 
			System.err.println(String.format("Ignored %d lines from %s", ignored, f.getName()));
		}
		
		br.close();
		return values;
	}
	
	private void writeLine(T modelObject, Vector<String> fieldOrder, PrintStream ps) {
		Class cls = modelObject.getClass();
		int fi = 0;
		for(String fieldName : fieldOrder) { 
			try {
				Field field = cls.getField(fieldName);
				Object value = field.get(modelObject);

				if(fi != 0) { ps.print("\t"); }
				if(value != null) { 
					ps.print(value.toString());
				} else { 
					ps.print("NA");
				}
				
			} catch (NoSuchFieldException e) {
				e.printStackTrace();
				ps.print("NA");
			} catch (IllegalAccessException e) {
				e.printStackTrace();
				ps.print("NA");
			}
			fi += 1;
		}
		if(fi > 0) { ps.println(); }
	}

	private static Pattern quotePattern = Pattern.compile("^\\s*\"(.*)\"\\s*$");

	private String extractQuoted(String quoted) { 
		Matcher m = quotePattern.matcher(quoted);
		if(m.matches()) { 
			String value = m.group(1);
			return value;
		} else { 
			return quoted;
		}
	}
	
	public T parseLine(String[] array, Vector<String> fieldOrder, Vector<Boolean> quoted) {
		
		if(fieldOrder.size() > array.length) {
			String arraystr = "";
			for(int i = 0; i < array.length; i++) { 
				arraystr += array[i] + " ";
			}
			String msg = String.format("fieldOrder.size() == %d (%s) exceeded array.length == %d : %s",
					fieldOrder.size(), fieldOrder.toString(), array.length, arraystr);
			throw new IllegalArgumentException(msg); 
		}
		
		T val = null;
		
		try {
			val = cls.newInstance();
			
			for(int i = 0; i < fieldOrder.size(); i++) { 
				String fieldName = fieldOrder.get(fieldOrder.size()-1-i);
				String valueString = array[array.length-1-i];
				if(quoted.get(i)) { 
					valueString = extractQuoted(valueString); 
				}
				
				boolean isNA = valueString.equals("NA");
			
				try { 
					Field f = cls.getField(fieldName);
					Class type = f.getType();

					if(Model.isSubclass(type, Double.class)) {
						try { 
							Double fieldValue = isNA ? null : Double.parseDouble(valueString);
							f.set(val, fieldValue);
						} catch(NumberFormatException nfe) { 
							f.set(val, null);  // missing value.
						}

					} else if(Model.isSubclass(type, Boolean.class)) { 
						Boolean fieldValue = isNA ? null : Boolean.parseBoolean(valueString);
						f.set(val, fieldValue);

					} else if (Model.isSubclass(type, Integer.class)) {
						try { 
							Integer fieldValue = isNA ? null : Integer.parseInt(valueString);
							f.set(val, fieldValue);
						} catch(NumberFormatException nfe) { 
							f.set(val, null);  // missing value.
						}

					} else if (Model.isSubclass(type, String.class)) { 
						String fieldValue = valueString;
						f.set(val, fieldValue);
					} else { 
						System.err.println(String.format(
								"Field %s has unsupported parsing type %s",
								f.getName(), type.getName()));
					}

				} catch (NoSuchFieldException e) {
				}
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		
		return val;
	}
}

class Test extends Model { 
	public String Weekend, Decision, Weather, Parents, Money;
}
