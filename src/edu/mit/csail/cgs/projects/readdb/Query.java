package edu.mit.csail.cgs.projects.readdb;

import org.apache.commons.cli.*;
import java.util.*;
import java.io.*;

/** 
 * Query hits from the database.  Reads chrom:start-stop:strand values from
 * stdin.  Repeats them on stdout along with the hit positions
 *
 * Usage:
 * Query --hostname nanog.csail.mit.edu --port 52000 --user foo --passwd bar [--quiet]
 *
 * --quiet means don't print any output.  This is useful for testing the query performance
 * without worrying about the time it takes to print the output.
 */

public class Query {

    private String alignname;
    private String hostname;
    private String username, password;
    private int portnum, histogram = -1;
    private boolean quiet, weights;
    

    public static void main(String args[]) throws Exception {
        Query query = new Query();
        query.parseArgs(args);
        query.run(System.in);
    }

    public void parseArgs(String args[]) throws IllegalArgumentException, ParseException {
        Options options = new Options();
        options.addOption("h","hostname",true,"server to connect to");
        options.addOption("P","port",true,"port to connect to");
        options.addOption("a","align",true,"alignment name");
        options.addOption("u","user",true,"username");
        options.addOption("p","passwd",true,"password");
        options.addOption("q","quiet",false,"quiet: don't print output");
        options.addOption("w","weights",false,"get and print weights in addition to positions");
        options.addOption("H","histogram",true,"produce a histogram with this binsize instead of printing all read positions");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        if (line.hasOption("port")) {
            portnum = Integer.parseInt(line.getOptionValue("port"));
        } else {
            portnum = -1;
        }
        if (line.hasOption("hostname")) {
            hostname = line.getOptionValue("hostname");
        } else {
            hostname = null;
        }
        if (line.hasOption("align")) {
            alignname = line.getOptionValue("align");
        } else {
            alignname = null;
        }
        if (line.hasOption("user")) {
            username = line.getOptionValue("user");
        } else {
            username = null;
        }
        if (line.hasOption("passwd")) {
            password = line.getOptionValue("passwd");
        } else {
            password = null;
        }
        if (line.hasOption("histogram")) {
            histogram = Integer.parseInt(line.getOptionValue("histogram"));
        }


        quiet = line.hasOption("quiet");
        weights = line.hasOption("weights");
    }

    public void run(InputStream instream) throws IOException, ClientException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
        String line;
        Client client;
        if (hostname != null && portnum != -1 && username != null && password != null) {
            client = new Client(hostname,
                                portnum,
                                username,
                                password);
        } else {
            client = new Client();
        }

        while ((line = reader.readLine()) != null) {
            try {
                String pieces[] = line.split("[\\:]");
                if (pieces.length != 3) {
                    System.err.println("Invalid query " + line);
                    continue;
                }
                String chr = pieces[0].replaceFirst("^chr","") + pieces[2];
                pieces = pieces[1].split("\\-");
                int start = Integer.parseInt(pieces[0]);
                int stop = Integer.parseInt(pieces[1]);
                if (histogram > 0) {
                    TreeMap<Integer,Integer> hits = client.getHistogram(alignname,
                                                                        chr,
                                                                        start,
                                                                        stop,
                                                                        histogram);
                    TreeMap<Integer,Float> weightsmap = null;
                    if (weights) {
                        weightsmap = client.getWeightHistogram(alignname,
                                                            chr,
                                                            start,
                                                            stop,
                                                            histogram);
                    }
                    if (!quiet) {
                        for (int i : hits.keySet()) {
                            if (weights) {
                                System.out.println(String.format("%d\t%d\t%f", i, hits.get(i), weightsmap.get(i)));  
                            } else {
                                System.out.println(String.format("%d\t%d", i, hits.get(i)));
                            }
                        }
                    }
                } else {
                    int[] hits = client.getHitsRange(alignname, chr, start,stop);
                    float[] w = null;
                    if (weights) {
                        w = client.getWeightsRange(alignname, chr, start,stop);
                    }
                    if (!quiet) {
                        System.out.println(line);
                        if (w != null) {
                            for (int i = 0; i < hits.length; i++) {
                                System.out.println(String.format("%d\t%.4f",hits[i], w[i]));
                            }
                        } else {
                            for (int i = 0; i < hits.length; i++) {
                                System.out.println(hits[i]);
                            }
                        }
                    }
                }
            } catch (Exception e) {
                System.err.println("FROM " + line);
                e.printStackTrace();
            }
        }

        client.close();
    }


    
}