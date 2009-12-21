package edu.mit.csail.cgs.tools.chipchip;

import java.io.*;
import java.util.*;
import java.sql.*;
import java.util.regex.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;

/**
 * MSP files are the output of the Young Lab's error model (also known
 * as the Rosetta error model.  
 *
 * Usage:
 * java edu.mit.csail.cgs.tools.chipchip.AddMSP --species "$SC;SGDv1" --analysis "Sc Gcn4 in YPD;linefitting normalization"  -- file1.msp file2.msp
 *
 */
public class AddMSP {

    private Organism species;
    private Genome genome;
    private String analysisname, analysisversion;
    private int analysisid;
    private java.sql.Connection core, chipchip;
    private int chromcol, poscol, ratiocol, xcol, pvalcol, pval3col, redcol, greencol, mediancol;
    private Map<String,Integer> chrommap;
    private ArrayList<String> fnames;

    public AddMSP() throws UnknownRoleException, SQLException {
        core = DatabaseFactory.getConnection("core");
        chipchip = DatabaseFactory.getConnection("chipchip");
        fnames = new ArrayList<String>();
    }

    public void parseArgs(String[] args) throws NotFoundException {
        species = null;
        genome = null;
        analysisname = null;
        analysisversion = null;
        int i;
        for (i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                String s = args[++i];
                String pieces[] = s.split(";");                
                species = new Organism(pieces[0]);
                if (pieces.length >= 2) {
                    genome = species.getGenome(pieces[1]);
                }
            }
            if (args[i].equals("--analysis") ||
                args[i].equals("--analysisname")) {
                String s = args[++i];
                String pieces[] = s.split(";");                
                analysisname = pieces[0];
                if (pieces.length >= 2) {
                    analysisversion = pieces[1];
                }
            }
            if (args[i].equals("--analysisversion")) {
                analysisversion = args[++i];
            }
            if (args[i].equals("--file")) {
                fnames.add(args[++i]);
            }
            if (args[i].equals("--")) {
                break;
            }
        }
        for (;i < args.length; i++) {
            fnames.add(args[i]);
        }

        for (i = 0; i < args.length; i++) {
            if (args[i].equals("--genome")) {
                genome = species.getGenome(args[++i]);
            }
        }
        if (species == null || genome == null) {
            throw new RuntimeException("Must supply --species 'species;genomeversion'");
        }
        if (analysisname == null || analysisversion == null) {
            throw new RuntimeException("Must supply --analysis 'analysisname;analysisversion'");
        }
    }

    /* parseArgs must be called first.  Ensures that the analysis specified
       by analysisname and analysisversion exists and sets analysisid accordingly. */
    public void createAnalysis() throws SQLException {
        Statement stmt = chipchip.createStatement();
        ResultSet rs = stmt.executeQuery("select id from rosettaanalysis where name = '" + 
                                         analysisname + "' and version = '" +
                                         analysisversion + "' and species =" + species.getDBID());
        if (rs.next()) {
            analysisid = rs.getInt(1);
        } else {
            rs.close();
            ArrayList<String> fieldnames = new ArrayList<String>();
            ArrayList<String> values = new ArrayList<String>();
            fieldnames.add("id"); values.add(Sequence.getInsertSQL(chipchip,"analysis_id"));
            fieldnames.add("name"); values.add("'" + analysisname + "'");
            fieldnames.add("version"); values.add("'" + analysisversion + "'");
            fieldnames.add("species"); values.add(Integer.toString(species.getDBID()));
            fieldnames.add("active"); values.add("1");
            Pattern exptname = Pattern.compile("\\s (\\w.*):(.*):(.*) vs (.*):(.*):(.*\\w)");
            Matcher matcher = exptname.matcher(analysisname);            
            if (matcher.matches()) {
                MetadataLoader loader = new MetadataLoader();
                fieldnames.add("factorone");
                values.add(Integer.toString(loader.getFactor(matcher.group(1)).getDBID()));
                fieldnames.add("cellsone");
                values.add(Integer.toString(loader.getCells(matcher.group(2)).getDBID()));
                fieldnames.add("conditionone");
                values.add(Integer.toString(loader.getCondition(matcher.group(3)).getDBID()));

                fieldnames.add("factortwo");
                values.add(Integer.toString(loader.getFactor(matcher.group(4)).getDBID()));
                fieldnames.add("cellstwo");
                values.add(Integer.toString(loader.getCells(matcher.group(5)).getDBID()));
                fieldnames.add("conditiontwo");
                values.add(Integer.toString(loader.getCondition(matcher.group(6)).getDBID()));               
            }

            StringBuffer cmd = new StringBuffer("insert into rosettaanalysis(");
            for (int i = 0; i < fieldnames.size(); i++) {
                cmd.append(fieldnames.get(i));
                if (i < fieldnames.size() -1 ) {
                    cmd.append(", ");
                }
            }
            cmd.append(") values(");
            for (int i = 0; i < values.size(); i++) {
                cmd.append(values.get(i));
                if (i < values.size() -1 ) {
                    cmd.append(", ");
                }
            }
            cmd.append(")");
            

            stmt.executeUpdate(cmd.toString());
            rs = stmt.executeQuery(Sequence.getLastSQLStatement(chipchip,"analysis_id"));
            rs.next();
            analysisid = rs.getInt(1);            
        }
        rs.close();

        rs = stmt.executeQuery("select count(*) from rosettaToGenome where analysis = " + 
                            analysisid + " and genome = " + genome.getDBID());
        rs.next();
        if (rs.getInt(1) == 0) {
            stmt.executeUpdate("insert into rosettaToGenome(analysis, genome) values (" +
                            analysisid + "," + genome.getDBID() + ")");
        }
        rs.close();
        stmt.close();
    }

    public void createChromMap() {
        chrommap = genome.getChromIDMap();
        if (species.getName().equals("Homo sapies")) {
            chrommap.put("23",chrommap.get("X"));
            chrommap.put("24",chrommap.get("Y"));
            chrommap.put("25",chrommap.get("mt"));
        }
        if (species.getName().equals("Mus musculus")) {
            chrommap.put("20",chrommap.get("X"));
            chrommap.put("21",chrommap.get("Y"));
            chrommap.put("22",chrommap.get("mt"));
        }
        if (species.getName().equals("Danio rerio")) {
            chrommap.put("26",chrommap.get("Un"));
            chrommap.put("27",chrommap.get("NA"));
        }
        if (species.getName().equals("Drosophila melanogaster")) {
            chrommap.put("1",chrommap.get("2L"));
            chrommap.put("2",chrommap.get("2R"));
            chrommap.put("3",chrommap.get("3L"));
            chrommap.put("4",chrommap.get("4"));
            chrommap.put("5",chrommap.get("3R"));
            chrommap.put("6",chrommap.get("X"));
            chrommap.put("7",chrommap.get("Y"));
        }

    }

    public void getColumnHeaders(String hline) {
        String colnames[] = hline.split("\\t");
        for (int i = 0; i < colnames.length; i++) {
            colnames[i] = colnames[i].toLowerCase();
            if (colnames[i].equals("chr")) {
                chromcol = i;
            } else if (colnames[i].equals("pos")) {
                poscol = i;
            } else if (colnames[i].equals("ratio")) {
                ratiocol = i;
            } else if (colnames[i].equals("x'")) {
                xcol = i;
            } else if (colnames[i].equals("pval1")) {
                pvalcol = i;
            } else if (colnames[i].equals("pval3")) {
                pval3col = i;
            } else if (colnames[i].equals("red")) {
                redcol = i;
            } else if (colnames[i].equals("green")) {
                greencol = i;
            } else if (colnames[i].equals("medianofratios")) {
                mediancol = i;
            }
        }       
    }

    public void readFile(String fname) throws IOException, SQLException {
        chipchip.setAutoCommit(false);
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        PreparedStatement insert = chipchip.prepareStatement("insert into rosettaresults(analysis, chromosome, position, ratio, " +
                                                             "x, pval, pval3, red, green, medianofratios) values (" +
                                                             "?,?,?,?,?,?,?,?,?,?)");
        insert.setInt(1,analysisid);
        String headerline = reader.readLine();
        getColumnHeaders(headerline);
        String line;     
        int added = 0;
        int duplicate = 0;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");
            if (pieces[mediancol].equals("FLAG") ||
                pieces[ratiocol].equals("FLAG") ||
                pieces[xcol].equals("FLAG")) {
                continue;
            }

            insert.setInt(2, chrommap.get(pieces[chromcol]));
            insert.setInt(3, Integer.parseInt(pieces[poscol]));
            insert.setDouble(4, Double.parseDouble(pieces[ratiocol]));
            insert.setDouble(5, Double.parseDouble(pieces[xcol]));
            insert.setDouble(6, Double.parseDouble(pieces[pvalcol]));
            insert.setDouble(7, Double.parseDouble(pieces[pval3col]));
            insert.setDouble(8, Double.parseDouble(pieces[redcol]));
            insert.setDouble(9, Double.parseDouble(pieces[greencol]));
            insert.setDouble(10, Double.parseDouble(pieces[mediancol]));
            try {
                insert.execute();
            } catch (SQLException e) {
                if (e.toString().matches(".*unique.*")) {
                    duplicate++;
                    continue;
                } else {
                    throw e;
                }
            }
            added++;
        }
        chipchip.commit();        
        System.err.println("Added " + added + " values from " + fname);
        System.err.println("Dropped " + duplicate + " apparently duplicate probe positions");
    }

    public static void main(String args[]) throws Exception {
        AddMSP msp = new AddMSP();
        msp.parseArgs(args);
        msp.createAnalysis();
        msp.createChromMap();
        for (String fname : msp.fnames) {
            msp.readFile(fname);
        }
    }
}
