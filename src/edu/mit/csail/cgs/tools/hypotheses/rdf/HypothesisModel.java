package edu.mit.csail.cgs.tools.hypotheses.rdf;

import java.sql.SQLException;
import java.util.*;

import com.hp.hpl.jena.db.DBConnection;
import com.hp.hpl.jena.db.IDBConnection;
import com.hp.hpl.jena.query.*;
import com.hp.hpl.jena.rdf.model.*;
import com.hp.hpl.jena.vocabulary.DC;
import com.hp.hpl.jena.vocabulary.VCARD;

public class HypothesisModel { 
    
    public static void main(String[] args) { 
        try {
            String className = "com.mysql.jdbc.Driver";
            Class.forName(className);
            
            String dbName = "jena";
            System.out.println("Connection to JENA database \"" + dbName + "\"");
            String url = "jdbc:mysql://nanog.csail.mit.edu/" + dbName;
            
            IDBConnection cxn = new DBConnection(url, "tdanford", "mysql4twd", "MySQL");
            
            ModelMaker maker = ModelFactory.createModelRDBMaker(cxn);
            Model model = maker.createDefaultModel();
            
            Resource rec = null;

            
            String cgsUrl = "http://cgs.csail.mit.edu/rdf";
            
            Property name = model.createProperty(cgsUrl, "#name");
            Property buildOf = model.createProperty(cgsUrl, "#buildof");
            
            model.createResource(cgsUrl + "/organism#mouse").
                addProperty(name, "Mus musculus");
            
            model.createResource(cgsUrl + "/genome#mm8").
                addProperty(name, "mm8").
                addProperty(buildOf, model.createResource(cgsUrl + "/genome#mouse"));
            
            String query = 
            	"select ?x " +
            	"where { ?x <" + name.getURI() + "> \"mm8\" } ";
            Query q = QueryFactory.create(query);
            QueryExecution exec = QueryExecutionFactory.create(q, model);
            
            ResultSet rs = exec.execSelect();
            while(rs.hasNext()) { 
            	QuerySolution sol = rs.nextSolution();
            	System.out.println("?x = " + sol.get("x"));
            }
            
            
            cxn.close();
            System.out.println("Finished.");
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
}