/*
 * Created on Apr 5, 2007
 */
package edu.mit.csail.cgs.datasets.expression.parsing;

import edu.mit.csail.cgs.datasets.expression.*;

import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;
import java.io.*;

public class TextTablePlatformParser {
    
    public static void main(String[] args) { 
        try {
			if(args.length != 3) { 
				System.err.println("USAGE: TextPlatformParser <platform_name> <platform_type> " + 
						"<platform_file>");
				System.err.println("\twhere <platform_type> is:");
				System.err.println("\t\t-1 : \"unknown type\"\n\t\t0: \"probe type\"\n" + 
						"\t\t1 : \"gene type\"");
				System.err.println("\tand <platform_file> has the format:\n" + 
						"\t\tcol 0: <probe_name>");
				System.exit(1);
			}

            TextTablePlatformParser parser = 
                new TextTablePlatformParser(args[0], Integer.parseInt(args[1]), new File(args[2]));
            System.out.println("Loaded " + parser.size() + " probes.");
            
            ExpressionInserter inserter = new ExpressionInserter();
            parser.insertIntoDB(inserter);
            
            inserter.close();
            System.out.println("Inserted probe-platform \"" + args[0] + "\" (" + args[1] + ")");
            
        } catch (IOException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    private Set<ProbeIdentifier> probes;
    private String platformName;
    private int platformType;
    
    public TextTablePlatformParser(String n, int t, File locFile) throws IOException {
        platformName = n;
        platformType = t;
        probes = new HashSet<ProbeIdentifier>();
        BufferedReader br = new BufferedReader(new FileReader(locFile));
        String line;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                ProbeIdentifier loc = new ProbeIdentifier(line);
                if(probes.contains(loc)) { throw new IllegalArgumentException("Duplicate probe: " + loc.name); }
                probes.add(loc);
            }
        }
        br.close();
    }
    
    public int size() { return probes.size(); }
    
    public void insertIntoDB(ExpressionInserter inserter) throws SQLException {
        TreeSet<String> probeNames = new TreeSet<String>();
        for(ProbeIdentifier ident : probes) { probeNames.add(ident.name); }

        inserter.beginTransaction();
        int platform = inserter.insertProbePlatform(platformName, platformType);
        inserter.insertProbes(probeNames, platform);
        inserter.commitTransaction();
    }
    
    private static class ProbeIdentifier {
        
        public String name;
        
        public ProbeIdentifier(String line) { 
            String[] array = line.split("\\s+");
            name = array[0];
        }
        
        public int hashCode() { 
            int code = 17;
            code += name.hashCode(); code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof ProbeIdentifier)) { return false; }
            ProbeIdentifier loc = (ProbeIdentifier)o;
            return loc.name.equals(name);
        }
    }
}
