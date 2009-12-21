package edu.mit.csail.cgs.datasets.fragdist;

import java.io.*;
import java.util.Map;
import java.util.HashMap;

public class CSV {
    
    public static Map<Double,Double> readFile(String fname) throws IOException, FileNotFoundException {
        return readFile(new BufferedReader(new FileReader(fname)));
    }
    public static Map<Double,Double> readFile(BufferedReader reader) throws IOException {
        HashMap<Double,Double> times = new HashMap<Double,Double>();
        boolean inDatapoints = false;
        String line;
        Double time = -1.0, intensity = -1.0;
        while ((line = reader.readLine()) != null) {
            if (line.matches(".*,.*")) {
                String pieces[] = line.split(",");
                if (pieces.length != 2) {
                    continue;
                }
                if (pieces[0].equals("Time") &&
                    pieces[1].equals("Value")) {
                    inDatapoints = true;
                } else if (inDatapoints) {
                    time = Double.parseDouble(pieces[0]);
                    intensity = Double.parseDouble(pieces[1]);
                    //                    System.err.println("PUTTING " + intensity + " for time " + time);
                    times.put(time,intensity);
                }                
            } else if (inDatapoints) {
                inDatapoints = false;
            }

        }
        return times;
    }

}
