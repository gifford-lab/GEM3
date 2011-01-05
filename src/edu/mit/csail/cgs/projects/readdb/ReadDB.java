package edu.mit.csail.cgs.projects.readdb;

import org.apache.commons.cli.*;
import java.util.*;
import java.io.*;

/**
 * <p>Generic command line client for the readdb
 * Usage:
 * <pre>ReadDB -h nanog.csail.mit.edu -P 5200 -u arolfe -p SECRET <command> .....</apre>
 *   omitting connection info uses default connection info in ~/.readdb_passwd
 *
 * <p>command can be:
 * <ul><li> exists alignname
 *  <li>getchroms alignname
 *  <li>getacl alignname
 *  <li>setacl alignname arolfe add write (add or delete; read, write or admin)
 *  <li>getcount alignname
 *  <li>getcount alignname chromname (eg, chromname = 1+)
 *  <li>addtogroup username groupname
 * 
 * <p>The --paired flag can be provided to make getweight, getcount, and getchroms work on paire-end rather than
 * single-end alignments
 *
 */
public class ReadDB {

    private Client client;
    private String[] otherargs;
    private boolean paired, isleft, noclose;

    public static void main(String args[]) {
        ReadDB readdb = null;
        try {
            readdb = new ReadDB();
            readdb.parseArgs(args);
            readdb.run();        
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (readdb != null && readdb.client != null && !readdb.noclose) {
                readdb.client.close();
            }
        }
    }

    public ReadDB() {
    }

    public void parseArgs(String args[]) throws IllegalArgumentException, ParseException, ClientException, IOException {
        Options options = new Options();
        options.addOption("H","hostname",true,"server to connect to");
        options.addOption("P","port",true,"port to connect to");
        options.addOption("u","user",true,"username");
        options.addOption("p","passwd",true,"password");
        options.addOption("d","paired",false,"work on paired alignment?");
        options.addOption("r","right",false,"query right side reads when querying paired alignments");
        options.addOption("C","noclose",false,"don't close the connection.  For debugging only");
        options.addOption("h","help",false,"print help");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        String hostname = null, username = null, password = null;
        int portnum = -1;
        if (line.hasOption("help")) {
            printHelp();
            System.exit(1);
        }
        if (line.hasOption("port")) {
            portnum = Integer.parseInt(line.getOptionValue("port"));
        }
        if (line.hasOption("hostname")) {
            hostname = line.getOptionValue("hostname");
        }
        if (line.hasOption("user")) {
            username = line.getOptionValue("user");
        }
        if (line.hasOption("passwd")) {
            password = line.getOptionValue("passwd");
        }
        if (portnum == -1 || hostname == null || username == null || password == null) {
            client = new Client();
        } else {
            client = new Client(hostname, portnum, username, password);        
        }
        paired = line.hasOption("paired");
        isleft = !line.hasOption("right");
        noclose = line.hasOption("noclose");
        otherargs = line.getArgs();
    }
    public void printHelp() {
        System.out.println("Administrative client for ReadDB");
        System.out.println("usage: java edu.mit.csail.cgs.projects.readdb.ReadDB command args...");
        System.out.println("  exists alignname");
        System.out.println("  getchroms alignname");
        System.out.println("  getacl alignname");
        System.out.println("  getcount alignname");
        System.out.println("  getcount alignname chromnameStrand   (eg, 1+)");
        System.out.println("  setacl alignname username|groupname add|delete write|read|admin ");
        System.out.println("  addtogroup username groupname");
    }

    public void run() throws IOException, ClientException {
        String cmd = otherargs[0];
        if (cmd.equals("shutdown")) {
            client.shutdown();
        } else if (cmd.equals("addtogroup")) {
            // username, groupname
            client.addToGroup(otherargs[1], otherargs[2]);

        } else {
            String align = otherargs[1];
            if (cmd.equals("exists")) {
                if (client.exists(align)) {
                    System.out.println("Exists");
                } else {
                    System.out.println("Doesn't exist");
                }
            } else {
                if (!client.exists(align)) {
                    System.err.println("No such alignment " + align);
                    return;
                }
                if (cmd.equals("getchroms")) {
                    for (Integer i : client.getChroms(align,paired,isleft)) {
                        System.out.println(i);
                    }
                } else if (cmd.equals("getacl")) {
                    Map<String,Set<String>> acls = client.getACL(align);
                    for (String acl : acls.keySet()) {
                        System.out.println(acl);
                        for (String p : acls.get(acl)) {
                            System.out.println("\t" + p);
                        }
                    }           
                } else if (cmd.equals("setacl")) {
                    String princ = otherargs[2];
                    String op = otherargs[3];
                    String acl = otherargs[4];                            
                    Set<ACLChangeEntry> changes = new HashSet<ACLChangeEntry>();
                    changes.add(new ACLChangeEntry(ACLChangeEntry.opCode(op),
                                                   ACLChangeEntry.aclCode(acl),
                                                   princ));
                    client.setACL(align,changes);
                } else if (cmd.equals("getcount")) {
                    int count = 0;
                    if (otherargs.length == 3) {
                        count = client.getCount(align,Integer.parseInt(otherargs[2]),paired,null,null,null,isleft,null);
                    } else {
                        count = client.getCount(align,paired,isleft,null);
                    }
                    if (otherargs.length == 3) {
                        System.err.println("Count in " + align + " chrom " + otherargs[2] + " is " + count);
                    } else {
                        System.err.println("Total count in " + align + " is " + count);
                    }
                } else if (cmd.equals("getweight")) {
                    double weight = 0;
                    if (otherargs.length == 3) {
                        weight = client.getWeight(align,Integer.parseInt(otherargs[2]),paired,null,null,null,isleft,null);
                    } else {
                        weight = client.getWeight(align,paired,isleft,null);
                    }
                    if (otherargs.length == 3) {
                        System.err.println("Weight in " + align + " chrom " + otherargs[2] + " is " + weight);
                    } else {
                        System.err.println("Total weight in " + align + " is " + weight);
                    }
                } else if (cmd.equals("deletesinglealign")) {
                    client.deleteAlignment(align,false);
                } else if (cmd.equals("deletepairedalign")) {
                    client.deleteAlignment(align,true);
                } else {
                    System.err.println("Unknown command " + cmd);
                }


            }
        }
    }
}