package edu.mit.csail.cgs.tools.genenames;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

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

public class GeneNames {
	
	public static PropertiesLoader props;
    public static String dbName;
    public static String host;
    public static String dbURL;
    public static String dbType;
    public static String cgsURL;
    public static String gnURL;    
    
    static {
    	props = new PropertiesLoader("geneprops");
    	dbName = props.getValue("dbName");
    	host = props.getValue("host");
    	dbURL = props.getValue("dbURL");
    	dbType = props.getValue("dbType");
    	cgsURL = props.getValue("cgsURL");
    	gnURL = props.getValue("gnURL");
    	
        try {
        	String driverName = props.getValue("driver");
            Class.forName(driverName);
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    private IDBConnection cxn;
    private Model model;
    
    private Property synonymOf, value, annotationOf;
    
    public GeneNames(String u, String p, boolean clear) { 
        cxn = new DBConnection(host+dbName, u, p, dbType);
        ModelMaker maker = ModelFactory.createModelRDBMaker(cxn);
        model = maker.createModel(dbURL);
        System.out.println("Connected to Model: " + dbURL);

        synonymOf = model.getProperty(gnURL + "/properties#synonymOf");
        value = model.getProperty(gnURL + "/properties#value");
        annotationOf = model.getProperty(gnURL + "/properties#annotationOf");
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

        synonymOf = model.getProperty(gnURL + "/properties#synonymOf");
        value = model.getProperty(gnURL + "/properties#value");
        annotationOf = model.getProperty(gnURL + "/properties#annotationOf");
    }
    
    public void addGeneSet(Genome g, String tablename) { 
        RefGeneGenerator<NamedRegion> gen = 
            new RefGeneGenerator<NamedRegion>(g, tablename);
        ChromRegionIterator chroms = new ChromRegionIterator(g);
        Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
        
        String geneURL = gnURL + "/genes#";
        int i = 0;

        while(genes.hasNext()) { 
            Gene gene = genes.next();
            String id = gene.getID();
            Resource idrec = model.createResource(geneURL + id);
            idrec.addProperty(value, id);
            
            for(String alias : gene.getAllNames()) { 
                if(!(alias.equals(id))) { 
                    Resource aliasrec = model.createResource(geneURL + alias);
                    aliasrec.addProperty(value, alias);
                    aliasrec.addProperty(synonymOf, idrec);
                    idrec.addProperty(synonymOf, aliasrec);
                }
            }
            
            i++;
            if(i % 100 == 0) { System.out.print("."); System.out.flush(); }
            if(i % 1000 == 0) { System.out.print("*"); System.out.flush(); }
        }
        
        System.out.println("\nAdded " + i + " genes.");
    }
    
    public ResultSet queryModel(String query) { 
        Query q = QueryFactory.create(query);
        QueryExecution qe = QueryExecutionFactory.create(q, model);
        ResultSet rs = qe.execSelect();
        return rs;
    }
    
    public static void main(String[] args) {
    	String username = props.getValue("user");
    	String password = props.getValue("password");
        GeneNames gn = new GeneNames(username, password, false);
        
        for(int i = 0; i < args.length; i++) { 
            if(args[i].equals("clear")) {
                gn.clear();
                
            } else if(args[i].startsWith("*")) { 
                String values = args[i].substring(1, args[i].length());
                String[] array = values.split(",");
                String genome = array[0], tablename = array[1];

                try {
                    gn.addGeneSet(Organism.findGenome(genome), tablename);
                } catch (NotFoundException e) {
                    e.printStackTrace();
                }

            } else { 
                String id = args[i];
                String query = "select ?synname ?syn ?target where { " +
                        "?syn <" + gn.synonymOf.getURI() + "> ?target . " +
                        "?target <" + gn.value.getURI() + "> \"" + id + "\" . " +
                        "?syn <" + gn.value.getURI() + "> ?synname }";
                //String query = "select ?synname where { ?synname <" + gn.synonym.getURI() + "> ?y }";

                ResultSet rs = gn.queryModel(query);
                
                while(rs.hasNext()) { 
                    QuerySolution soln = rs.nextSolution();
                    System.out.println("Name: " + soln.get("synname"));
                    System.out.println("\t" + soln.get("syn") + " -> " + soln.get("target"));
                }                
            }
        }
        
        gn.close();
    }
}
