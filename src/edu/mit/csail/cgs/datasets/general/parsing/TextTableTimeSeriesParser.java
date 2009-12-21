/**
 * April 7th, 2007
 * 
 * @author Timothy Danford
 */
package edu.mit.csail.cgs.datasets.general.parsing;

import java.io.*;
import java.sql.SQLException;
import java.util.Vector;
import edu.mit.csail.cgs.datasets.general.*;

public class TextTableTimeSeriesParser { 

	public static void main(String[] args) { 
		try { 
			TextTableTimeSeriesParser parser = new TextTableTimeSeriesParser(new File(args[0]));
			TimeSeriesLoader loader = new TimeSeriesLoader();
			parser.insertIntoDB(loader);
			loader.close();
		} catch(SQLException e) { 
			e.printStackTrace(System.err);
		} catch(IOException e) { 
			e.printStackTrace(System.err);
		}
	}

	private String name;
	private Vector<String> points;

	public TextTableTimeSeriesParser(File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));

		String line;
		name = null;
		points = new Vector<String>();

		while((line = br.readLine()) != null) { 
			line = line.trim();
			if(line.length() > 0) { 
				if(name == null) { 
					name = line;
				} else { 
					points.add(line);
				}
			}
		}

		br.close();
	}

	public void insertIntoDB(TimeSeriesLoader loader) throws SQLException { 
		loader.beginTransaction();

		int tsid = loader.insertTimeSeries(name);
		TimeSeries series = loader.loadTimeSeries(tsid);
		for(int i = 0; i < points.size(); i++) { 
			int tpid = loader.insertTimePoint(series, points.get(i), i);
		}

		loader.commitTransaction();
		System.out.println("Loaded " + points.size() + " points in Time Series \"" + name + "\"");
	}
}
