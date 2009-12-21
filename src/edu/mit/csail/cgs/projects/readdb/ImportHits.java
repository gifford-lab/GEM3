package edu.mit.csail.cgs.projects.readdb;

import org.apache.commons.cli.*;
import java.net.*;
import java.util.*;
import java.io.*;

/**
 * Imports hits to the db.
 * Usage:
 * cat hits.txt | ImportHits -h nanog.csail.mit.edu -P 5200 -a "Gcn4 ChipSeq" -u arolfe -p SECRET
 *
 * Lines in the input must be of the form
 * chromosome\tstart\tstrand\tweight
 *
 * where start is the position of the 5' end of the read
 * 
 *
 * ImportHits strips "chr" from the front of the
 * chromosome name and appends " +" or " -"
 * to indicate the strandedness of the hits stored.  The position stored in the DB is
 * the average of start and stop.
 */
public class ImportHits {

    String alignname;
    String hostname;
    String username, password;
    int portnum;

    public static void main(String args[]) throws Exception {
        ImportHits importer = new ImportHits();
        importer.parseArgs(args);
        importer.run(System.in);
    }

    public ImportHits(String hostname,
                      int port,
                      String alignname,
                      String username, String password) {
        this.hostname = hostname;
        this.portnum = port;
        this.alignname = alignname;
        this.username = username;
        this.password = password;
    }

    /* use this constructor and then call parseArgs */
    public ImportHits () { 
        username = null; 
        password = null;
        hostname = null;
        portnum = -1;
    }    
    public void parseArgs(String args[]) throws IllegalArgumentException, ParseException {
        Options options = new Options();
        options.addOption("h","hostname",true,"server to connect to");
        options.addOption("P","port",true,"port to connect to");
        options.addOption("a","align",true,"alignment name");
        options.addOption("u","user",true,"username");
        options.addOption("p","passwd",true,"password");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        if (line.hasOption("port")) {
            portnum = Integer.parseInt(line.getOptionValue("port"));
        }
        if (line.hasOption("hostname")) {
            hostname = line.getOptionValue("hostname");
        }
        if (line.hasOption("align")) {
            alignname = line.getOptionValue("align");
        } else {
            System.err.println("Must supply alignment name as --align");
            throw new IllegalArgumentException("Must supply alignment name as --align");
        }        
        if (line.hasOption("user")) {
            username = line.getOptionValue("user");
        }
        if (line.hasOption("passwd")) {
            password = line.getOptionValue("passwd");
        }


    }

    public void run(InputStream instream) throws IOException, ClientException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
        Map<String,ArrayList<Integer>> posmap = new HashMap<String,ArrayList<Integer>>();
        Map<String,ArrayList<Float>> weightmap = new HashMap<String,ArrayList<Float>>();
        String line;
        int lineno = 0;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");            
            String chr = pieces[0].replaceFirst("^chr","") + pieces[2];
            int pos = Integer.parseInt(pieces[1]);
            float weight = Float.parseFloat(pieces[3]);


            if (!posmap.containsKey(chr)) {
                posmap.put(chr, new ArrayList<Integer>());
                weightmap.put(chr, new ArrayList<Float>());
            }
            posmap.get(chr).add(pos);
            weightmap.get(chr).add(weight);
            if (lineno++ % 100000 == 0) {
                System.err.println("Read through line " + lineno);
            }

        }
        Client client;
        if (hostname != null && portnum > 0 && username != null && password != null) {
            client = new Client(hostname,
                                portnum,
                                username,
                                password);   
        } else {
            client = new Client();
        }

        for (String chr : posmap.keySet()) {
            List<Integer> poslist = posmap.get(chr);
            List<Float> weightlist = weightmap.get(chr);
            int[] hits = new int[poslist.size()];
            float[] weights = new float[poslist.size()];
            for (int i = 0; i < hits.length; i++) {
                hits[i] = poslist.get(i);
                weights[i] = weightlist.get(i);
            }
            System.err.println("Storing " + chr + " numhits=" + hits.length);
            client.store(alignname, chr, hits,weights);
        }       
        client.close();
    }

}
