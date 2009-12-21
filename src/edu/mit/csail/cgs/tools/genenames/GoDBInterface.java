package edu.mit.csail.cgs.tools.genenames;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.util.*;

import com.hp.hpl.jena.db.DBConnection;
import com.hp.hpl.jena.db.IDBConnection;
import com.hp.hpl.jena.query.*;
import com.hp.hpl.jena.rdf.model.*;
import com.hp.hpl.jena.vocabulary.DC;
import com.hp.hpl.jena.vocabulary.VCARD;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;

public class GoDBInterface {
    
    public static PropertiesLoader props;
    public static String dbName;
    public static String host;
    public static String dbURL;
    public static String dbType;

    static { 
    	props = new PropertiesLoader("goprops");
    	dbName = props.getValue("dbName");
    	host = props.getValue("host");
    	dbURL = props.getValue("dbURL");
    	dbType = props.getValue("dbType");
    	
        try {
        	String driver = props.getValue("driver");
            Class.forName(driver);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    private IDBConnection cxn;
    private Model model;
    
    private Property synonymOf, value, annotationOf;
    
    public GoDBInterface(String u, String p, boolean clear) { 
        cxn = new DBConnection(host+dbName, u, p, dbType);
        ModelMaker maker = ModelFactory.createModelRDBMaker(cxn);
        model = maker.createModel(dbURL);
        System.out.println("Connected to Model: " + dbURL);
    }
    
    public void close() { 
        try {
            cxn.close();
            cxn = null;
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public void clear() { 
        model.removeAll();
    }

    public ResultSet queryModel(String query) { 
        Query q = QueryFactory.create(query);
        QueryExecution qe = QueryExecutionFactory.create(q, model);
        ResultSet rs = qe.execSelect();
        return rs;
    }
    
    public static void main(String[] args) { 
        GoDBInterface gn = new GoDBInterface("tdanford", "mysql4twd", false);
        
        if(args.length > 0) { 
            if(args[0].equals("clear")) {
                System.out.println("Clearing the model.");
                gn.clear();
            } else {
                gn.model.begin();
                //System.out.println("Clearing the model.");
                //gn.clear();

                String fname = args[0];
                File inputFile = new File(fname);
                try {
                    System.out.println("Reading: " + fname);
                    Model mm = ModelFactory.createDefaultModel();
                    InputStream ins = new FileInputStream(inputFile);
                    mm.read(ins, "");
                    ins.close();
                    System.out.println("Finished.");
                    
                    System.out.println("Inputting the Model...");
                    gn.model.add(mm);
                    System.out.println("Finished.");
                    
                    gn.model.commit();
                    
                } catch (IOException e) {
                    gn.model.abort();
                    e.printStackTrace();
                }
            }
        }
        
        InputStream is = System.in;
        BufferedReader br = new BufferedReader(new InputStreamReader(is));
        
        try { 
            String line;
            String cmd = "";
            while((line = br.readLine()) != null) { 
                line = line.trim();
                if(line.endsWith("\\\\")) { 
                    cmd += line + " ";
                } else { 
                    String query = cmd;
                    cmd = "";
                    
                    System.out.println("Query: \"" + query + "\"");
                    
                    ResultSet rs = gn.queryModel(query);
                    List vars = rs.getResultVars();
                    while(rs.hasNext()) { 
                        QuerySolution soln = rs.nextSolution();
                        Iterator varitr = vars.iterator();
                        while(varitr.hasNext()) { 
                            String var = (String)varitr.next();
                            RDFNode node = soln.get(var);
                            System.out.println(var + "=" + node);
                        }
                        System.out.println();
                    }                   
                }
            }
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
        
        gn.close();
    }
}

