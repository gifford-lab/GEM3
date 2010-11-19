package edu.mit.csail.cgs.projects.readdb;

import org.apache.commons.cli.*;
import java.util.*;
import java.io.*;

/**
 * Imports hits to the db.
 * Usage:
 * cat hits.txt | ImportHits -H nanog.csail.mit.edu -P 5200 -a "Gcn4 ChipSeq" -u arolfe -p SECRET
 *
 * Lines in the input must be of the form
 * chromosome\tstart\tstrand\tlength\tweight
 *  or
 * chromosomeone\tstartone\tstrandone\tlengthone\tchromosometwo\tstarttwo\tstrandtwo\tlengthtwo\tweight
 *
 * where start is the position of the 5' end of the read.
 *
 * For CSAIL use:
 * - The chromosome must be numeric and should be the chromosome id from core.  
 * - The alignment should be numeric and should the alignment identifier from the chipseq schema
 *
 */
public class ImportHits {

    String alignname;
    String hostname;
    String username, password;
    int portnum;
    private Client client;
    private int chunk = 10000000;

    public static void main(String args[])  {
        ImportHits importer = null;
        try {
            importer = new ImportHits();
            importer.parseArgs(args);
            importer.run(System.in);
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (importer != null && importer.client != null) {
                importer.client.close();
            }
        }
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
        options.addOption("H","hostname",true,"server to connect to");
        options.addOption("P","port",true,"port to connect to");
        options.addOption("a","align",true,"alignment name");
        options.addOption("u","user",true,"username");
        options.addOption("p","passwd",true,"password");
        options.addOption("h","help",false,"print help message");
        options.addOption("c","chunk",true,"send this many hits to the server at once");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        if (line.hasOption("help")) {
            printHelp();
            System.exit(0);
        }

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
        if (line.hasOption("chunk")) {
            chunk = Integer.parseInt(line.getOptionValue("chunk"));
        }

    }
    public void printHelp() {
        System.out.println("ImportHits to ReadDB");
        System.out.println("usage: cat foo.sam | java edu.mit.csail.cgs.projects.readdb.SAMToReadDB | java edu.mit.csail.cgs.projects.readdb.ImportHits \\");
        System.out.println(" --align alignmentname");
        System.out.println(" [--help] print usage");
        System.out.println("");
        System.out.println("Input format is tab delimited with either five or nine fields per line.");
        System.out.println("For single-ended reads, fields are ");
        System.out.println("  (1) chromosome (2) position of 5' end of read (3) strand (4) length (5) weight");
        System.out.println("For paired-end reads, fields are ");
        System.out.println("  (1) L chromosome (2) position of 5' end of L read (3) L strand (4) L length ");
        System.out.println("  (5) R chromosome (6) position of 5' end of R read (7) R strand (8) R length (9) weight");
    }

    public void run(InputStream instream) throws IOException, ClientException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
        String line;
        int lineno = 0;
        List<SingleHit> hits = new ArrayList<SingleHit>();
        List<PairedHit> paired = new ArrayList<PairedHit>();
        if (hostname != null && portnum > 0 && username != null && password != null) {
            client = new Client(hostname,
                                portnum,
                                username,
                                password);   
        } else {
            client = new Client();
        }
        System.err.println("Created Client");
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");            
            if (pieces.length == 5) {
                hits.add(new SingleHit(Integer.parseInt(pieces[0]),
                                       Integer.parseInt(pieces[1]),
                                       Float.parseFloat(pieces[4]),
                                       pieces[2].equals("+"),
                                       Short.parseShort(pieces[3])));
            } else if (pieces.length == 9) {
                paired.add(new PairedHit(Integer.parseInt(pieces[0]),
                                         Integer.parseInt(pieces[1]),
                                         pieces[2].equals("+"),
                                         Short.parseShort(pieces[3]),
                                         Integer.parseInt(pieces[4]),
                                         Integer.parseInt(pieces[5]),
                                         pieces[6].equals("+"),
                                         Short.parseShort(pieces[7]),
                                         Float.parseFloat(pieces[8])));
            } else {
                System.err.println("Bad line size " + line);
            }
            if (lineno++ % 100000 == 0) {
                System.err.println("Read through line " + lineno);
            }
            if (lineno % chunk == 0) {
                if (hits.size() > 0) {
                    try {
                        client.storeSingle(alignname, hits);
                        hits.clear();
                    } catch (Exception e) {
                        System.err.println("Failed: " + e.toString());
                        e.printStackTrace();
                    }
                    
                }
                if (paired.size() > 0) {
                    try {
                        client.storePaired(alignname, paired);
                        paired.clear();
                    } catch (Exception e) {
                        System.err.println("Failed: " + e.toString());
                        e.printStackTrace();
                    }
                }
                
            }


        }
        System.err.println("Read lines");
        if (hits.size() > 0) {
            try {
                client.storeSingle(alignname, hits);
            } catch (Exception e) {
                System.err.println("Failed: " + e.toString());
                e.printStackTrace();
            }

        }
        if (paired.size() > 0) {
            try {
                client.storePaired(alignname, paired);
            } catch (Exception e) {
                System.err.println("Failed: " + e.toString());
                e.printStackTrace();
            }
        }
        System.err.println("Stored");
    }

}
