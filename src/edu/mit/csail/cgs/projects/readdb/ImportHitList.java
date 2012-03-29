package edu.mit.csail.cgs.projects.readdb;

import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class ImportHitList {

	String alignname;
	String hostname;
	String username, password;
	int portnum;
	private Client client;
	private int chunk = 10000000;

	public ImportHitList(String hostname,
			int port,
			String alignname,
			String username, String password) {
		this.hostname = hostname;
		this.portnum = port;
		this.alignname = alignname;
		this.username = username;
		this.password = password;
	}
	
	public ImportHitList() {
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
        options.addOption("c","chunk",true,"send this many hits to the server at once");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, true );            
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

	public void importHitList(List<PairedHit> paired) throws IOException, ClientException {
		if (hostname != null && portnum > 0 && username != null && password != null) {
			client = new Client(hostname,
					portnum,
					username,
					password);   
		} else {
			client = new Client();
		}
		System.err.println("Created Client");
		for (int i=0; i<paired.size()/chunk; i++) {
			try {
				client.storePaired(alignname, paired.subList(i*chunk,(i+1)*chunk));
				System.err.println("Stored "+((i+1)*chunk));
			} catch (Exception e) {
				System.err.println("Failed: " + e.toString());
				e.printStackTrace();
			}
		}
		try {
			client.storePaired(alignname, paired.subList((paired.size()/chunk)*chunk,paired.size()));
		} catch (Exception e) {
			System.err.println("Failed: " + e.toString());
			e.printStackTrace();
		}
		System.err.println("Stored all");
		client.close();
	}

}
