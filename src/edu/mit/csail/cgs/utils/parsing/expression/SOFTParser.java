package edu.mit.csail.cgs.utils.parsing.expression;

import java.util.*;
import java.util.regex.*;

import java.io.*;

import edu.mit.csail.cgs.utils.Pair;

public class SOFTParser {
	
	public static void main(String[] args) {
		File f = new File(args[0]);
		try {
			SOFTParser parser = new SOFTParser(f);
			parser.printSummary();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private SOFTDatabase database;
	private SOFTPlatform platform;
	private SOFTSeries series;
	private Vector<SOFTSample> samples;
	
	private SOFTSample currentSample;
	private SOFTPlatform currentPlatform;
	private SOFTAttributes currentAttrs;
	
	private File output;
	private PrintStream outps;
	
	public SOFTParser(File f) throws IOException { 
		this(f, null);
	}
	
	public SOFTParser(File f, File outf) throws IOException { 
		database = null;
		platform = null;
		series = null;
		currentSample = null;
		currentPlatform = null;
		currentAttrs = null;
		output = outf;
		outps = null;
		samples = new Vector<SOFTSample>();

		parseFile(f);
	}
	
	public void parseFile(File f) throws IOException {
		if(output != null) { 
			outps = new PrintStream(new FileOutputStream(output));
		} else { 
			outps = null;
		}
		String line = null;
		BufferedReader br = new BufferedReader(new FileReader(f));
		
		while((line = br.readLine()) != null) { 
			String toParse = null;
			
			if((toParse = parseCaretLine(line)) != null) {
				Pair<String,String> kv = parseKeyValue(toParse);
				
				if(kv == null) { 
					String msg = String.format("Couldn't parse key/value pair: %s",
							toParse);
					throw new IllegalArgumentException(msg);
				}
				String key = kv.getFirst().trim(), value = kv.getLast();
				
				if(key.equals("DATABASE")) { 
					database = new SOFTDatabase();
					currentAttrs = database.getAttributes();
					currentSample = null;
					currentPlatform = null;
					
					System.out.println(String.format("DATABASE: %s", value));

				} else if (key.equals("PLATFORM")) { 
					platform = new SOFTPlatform();
					currentAttrs = platform.getAttributes();
					currentSample = null;
					currentPlatform = platform;

					System.out.println(String.format("PLATFORM: %s", value));

				} else if (key.equals("SERIES")) {
					series = new SOFTSeries();
					currentAttrs = series.getAttributes();
					currentSample = null;
					currentPlatform = null;

					System.out.println(String.format("SERIES: %s", value));

				} else if (key.equals("SAMPLE")) {

					if(outps != null) { 
						if(currentPlatform != null) { 
							currentPlatform.printSampleHeaderLine(outps);
						} else if (currentSample != null) { 
							currentSample.printSampleLine(platform, outps);
						}
					}
					
					SOFTSample newSample = new SOFTSample();
					
					if(outps == null) { 
						samples.add(newSample);
					}
					
					currentAttrs = newSample.getAttributes();
					currentSample = newSample;
					currentPlatform = null;

					System.out.println(String.format("SAMPLE: %s", value));

				} else { 
					String msg = String.format("Unknown SOFT file type: \"%s\"", key);
					throw new IllegalArgumentException(msg);
				}

				currentAttrs.addKeyValue(key, value);

			} else if((toParse = parseBangLine(line)) != null) {
				
				if(toParse.equals("platform_table_begin")) {
					if(currentPlatform==null) {
						String msg = String.format("Unexpected platform_table_begin");
						throw new IllegalArgumentException(msg);
					}

				} else if(toParse.equals("platform_table_end")) { 
					if(currentPlatform==null) {
						String msg = String.format("Unexpected platform_table_end");
						throw new IllegalArgumentException(msg);
					}
					
					System.out.println(String.format("# Platform Rows: %d", 
							currentPlatform.getTableSize()));

				} else if(toParse.equals("sample_table_begin")) {
					if(currentSample==null) {
						String msg = String.format("Unexpected sample_table_begin");
						throw new IllegalArgumentException(msg);
					}

				} else if(toParse.equals("sample_table_end")) {
					if(currentSample==null) {
						String msg = String.format("Unexpected sample_table_end");
						throw new IllegalArgumentException(msg);
					}

					System.out.println(String.format("# Sample Rows: %d", 
							currentSample.getTableSize()));

				} else { 
					Pair<String,String> kv = parseKeyValue(toParse);
					if(kv == null) { 
						String msg = String.format("Couldn't parse key/value pair: %s",
								toParse);
						//throw new IllegalArgumentException(msg);
					} else { 
						String key = kv.getFirst(), value = kv.getLast();

						if(currentAttrs == null) { 
							String msg = String.format("Unexpected ! line: %s", line);
							throw new IllegalArgumentException(msg);
						}

						currentAttrs.addKeyValue(key, value);
					}
				}
				
			} else if((toParse = parseHashLine(line)) != null) {
				Pair<String,String> kv = parseKeyValue(toParse);
				if(kv == null) { 
					String msg = String.format("Couldn't parse key/value pair: %s",
							toParse);
					throw new IllegalArgumentException(msg);
				}
				String key = kv.getFirst(), value = kv.getLast();

				if(currentPlatform != null) {
					currentPlatform.getTableAttributes().addKeyValue(key, value);
					
				} else if (currentSample != null) { 
					currentSample.getTableAttributes().addKeyValue(key, value);
					
				} else { 
					String msg = String.format("Unexpected # line " +
							"not in platform or sample: \"%s\"", line);
					//throw new IllegalArgumentException(msg);
				}
				
			} else {
				if(currentPlatform != null) {
					if(currentPlatform.getHeader() == null) { 
						currentPlatform.setHeader(line);
					} else { 
						currentPlatform.addToTable(line);
					}
					
				} else if (currentSample != null) { 
					if(currentSample.getHeader() == null) { 
						currentSample.setHeader(line);
					} else { 
						currentSample.addToTable(line);
					}
					
				} else { 
					String msg = String.format("Unexpected non-escaped line: \"%s\"", line);
					//throw new IllegalArgumentException(msg);
				}
			}
		}
		
		if(currentSample != null && outps != null) { 
			currentSample.printSampleLine(platform,outps);
		}
		
		br.close();
		if(outps != null) { 
			outps.close();
			outps = null;
		}
	}
	
	public SOFTDatabase getDatabase() { return database; }
	public SOFTPlatform getPlatform() { return platform; }
	public SOFTSeries getSeries() { return series; }
	public int getNumSamples() { return samples.size(); }
	public SOFTSample getSample(int i) { return samples.get(i); }
	
	public void printSummary() { 
		if(database != null) { 
			System.out.println(String.format("DATABASE: %s", 
					database.getAttributes().getCompleteValue("DATABASE")));
		}
		if(platform!= null) { 
			System.out.println(String.format("PLATFORM: %s", 
					platform.getAttributes().getCompleteValue("PLATFORM")));
		}
		if(series!= null) { 
			System.out.println(String.format("SERIES: %s", 
					series.getAttributes().getCompleteValue("SERIES")));
		}
		
		for(SOFTSample sample : samples) { 
			System.out.println(String.format("SAMPLE: %s", 
					sample.getAttributes().getCompleteValue("SAMPLE")));			
		}
	}

	public static Pattern caretLine = Pattern.compile("^\\^(.+)");
	public static Pattern bangLine = Pattern.compile("^!(.+)");
	public static Pattern hashLine = Pattern.compile("^#(.+)");
	
	public static Pattern keyValuePattern = Pattern.compile("\\s*([^=]+)\\s*=\\s*(.*)");
	
	public static Pair<String,String> parseKeyValue(String line) { 
		Matcher m = keyValuePattern.matcher(line);
		if(m.matches()) { 
			return new Pair<String,String>(m.group(1), m.group(2));
		} else { 
			return null;
		}
	}
	
	public static String parseCaretLine(String line) { 
		Matcher m = caretLine.matcher(line);
		if(m.matches()) { 
			return m.group(1);
		} else { 
			return null;
		}
	}
	
	public static String parseBangLine(String line) { 
		Matcher m = bangLine.matcher(line);
		if(m.matches()) { 
			return m.group(1);
		} else { 
			return null;
		}
	}
	
	public static String parseHashLine(String line) { 
		Matcher m = hashLine.matcher(line);
		if(m.matches()) { 
			return m.group(1);
		} else { 
			return null;
		}
	}
}
