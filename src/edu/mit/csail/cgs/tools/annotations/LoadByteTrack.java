package edu.mit.csail.cgs.tools.annotations;

import java.sql.*;
import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;

/* reads a set of FASTA-like files that contain one number (space separated) representing 
   one byte per base.  Creates table in the appropriate UCSC database for each file and 
   then uploads the data.  

   This is primarily intended for data about assembly quality, which is reported as one byte per base of the assembly */


public class LoadByteTrack {

    private static Organism organism;
    private static  Genome genome;
    private static java.sql.Connection cxn;
    private static ArrayList<String> fnames;
    private static String tablename;
    private static PreparedStatement insert, append;

    public static void parseArgs(String args[]) throws NotFoundException, SQLException  {
        String org = null, gen = null;
        tablename = null;
        fnames = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                String pieces[] = args[++i].split(";");
                if (pieces.length > 0) {
                    org = pieces[0];
                }
                if (pieces.length > 1) {
                    gen = pieces[1];
                }
            } 
            if (args[i].equals("--genome")) {
                gen = args[++i];
            }
            if (args[i].equals("--tablename")) {
                tablename = args[++i];
            }

            if (args[i].equals("--")) {
                i++;
                fnames = new ArrayList<String>();
                while (i < args.length) {
                    fnames.add(args[i++]);
                }
            }           
        }
        if (tablename == null) {
            throw new IllegalArgumentException("must supply --tablename to tell me how to name the table");
        }
        if (fnames == null) {
            throw new IllegalArgumentException("must supply filenames after -- ");
        }
        organism = new Organism(org);
        genome = organism.getGenome(gen);
        cxn = genome.getUcscConnection();
    }

    public static void makeTable() throws SQLException {
        Statement stmt = cxn.createStatement();
        try {
            if (DatabaseFactory.isOracle(cxn)) {
                stmt.execute("create table " + tablename + "(chromosome varchar2(300), quality blob)");
            } else if (DatabaseFactory.isMySQL(cxn)) {
                stmt.execute("create table " + tablename + "(chromosome varchar(300), quality longblob)");
            } else {
                throw new RuntimeException("Don't understand database type " + DatabaseFactory.getType(cxn));
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        insert = cxn.prepareStatement("insert into " + tablename + " (chromosome, quality) values (?,?)");
        append = cxn.prepareStatement("update " + tablename + " set quality = concat(quality,?) where chromosome = ?");
    }

    public static void processFile(String fname) throws IOException, SQLException  {
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String header = reader.readLine();
        header = header.replaceAll("^>","");
        header = header.replaceAll(".f(ast)?a","");
        header = header.replaceAll(" quality scores","");
        System.err.println("Header is " + header);
        String line;
        ArrayList<String[]> lines = new ArrayList<String[]>();
        byte[] bytes = new byte[0];
        insert.setString(1,header);
        insert.setBytes(2,bytes);
        insert.execute();
        System.err.println("Inserting " + header);
        append.setString(2,header);
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\s+");
            bytes = new byte[pieces.length];
            for (int i = 0; i < pieces.length; i++) {                
                if (pieces[i].length() == 0) {
                    break;
                }
                bytes[i] = Byte.parseByte(pieces[i]);                
            }
            append.setBytes(1,bytes);
            append.execute();
        }
        reader.close();
    }

    public static void main(String args[]) throws Exception {
        parseArgs(args);
        makeTable();
        for (String fname : fnames) {
            processFile(fname);
        }

    }

}
