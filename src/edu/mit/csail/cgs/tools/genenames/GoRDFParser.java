/*
 * Created on Jul 23, 2007
 */
package edu.mit.csail.cgs.tools.genenames;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import com.hp.hpl.jena.db.DBConnection;
import com.hp.hpl.jena.db.IDBConnection;
import com.hp.hpl.jena.query.*;
import com.hp.hpl.jena.rdf.model.*;

import com.hp.hpl.jena.rdf.model.*;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class GoRDFParser {
    
    public static void main(String[] args) {
		String fname = args[0];
		String versionName = args[1];
        System.out.println("Loading: " + fname);
        
        Model model = ModelFactory.createDefaultModel();
        
        File f = new File(fname);
        InputStream is = null;
        try {
            is = new FileInputStream(f);
            
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(0);
        }

        model.read(is, "");

        try {
            is.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        System.out.println("Loaded model.");

        String gobase = "http://www.geneontology.org/dtds/go.dtd#";
        
        Property acc = model.createProperty(gobase + "accession");
        Property name = model.createProperty(gobase + "name");
        Property isa = model.createProperty(gobase + "is_a");
        
        String query = "select ?uri ?acc ?name where { " +
        "?uri <" + acc.getURI() + "> ?acc . " +
        "?uri <" + name.getURI() + "> ?name }";
        
        String query2 = "select ?childacc ?parentacc where { " + 
        	"?child <" + isa.getURI() + "> ?parent . " +
        	"?child <" + acc.getURI() + "> ?childacc . " +
        	"?parent <" + acc.getURI() + "> ?parentacc }";

        Query q = QueryFactory.create(query);
        QueryExecution qe = QueryExecutionFactory.create(q, model);
        ResultSet rs = qe.execSelect();
        
        try {
			DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
			loader.beginTransaction();
			
			loader.insertVersion(versionName);
			FunctionVersion version = loader.getVersion(versionName);
			
			Map<String,Category> cats = new HashMap<String,Category>();
			
	        while(rs.hasNext()) { 
	        	QuerySolution soln = rs.nextSolution();
	        	RDFNode urinode = soln.getResource("uri");
	        	Literal namelit = soln.getLiteral("name");
	        	Literal acclit = soln.getLiteral("acc");
	        	
	        	String cn = namelit.toString();
	        	String cacc = acclit.toString();
	        	
	        	loader.insertCategory(version, cacc, cn);
	        	cats.put(cacc, loader.getCategory(version, cacc));
	        }
	        
	        System.out.println("Inserted " + cats.size() + " categories.");
	        
	        Query q2 = QueryFactory.create(query2);
	        QueryExecution qe2 = QueryExecutionFactory.create(q2, model);
	        rs = qe2.execSelect();
	        
	        int cr = 0;
	        while(rs.hasNext()) { 
	        	QuerySolution qs = rs.nextSolution();
	        	Literal childacclit = qs.getLiteral("childacc");
	        	Literal parentacclit = qs.getLiteral("parentacc");
	        	String childacc = childacclit.toString();
	        	String parentacc = parentacclit.toString();
	        	Category child = cats.get(childacc), parent = cats.get(parentacc);
	        	
				if(!cats.containsKey(childacc)) { System.err.println("No acc: " + childacc); }
				if(!cats.containsKey(parentacc)) { System.err.println("No acc: " + parentacc); }
	        	
	        	if(child != null && parent != null) { 
	        		loader.insertCategoryRelationship(child, parent);
		        	cr += 1;
		        	if(cr % 100 == 0) { System.out.print("."); System.out.flush(); }
		        	if(cr % 1000 == 0) { System.out.print("*"); System.out.flush(); }
	        	}
	        }
	        System.out.println("\nInserted " + cr + " relationships");
	        
	        for(int i = 2; i < args.length; i++) { 
	        	Vector<String[]> arrays = loadPairs(new File(args[i]));
	        	for(String[] array : arrays) { 
	        		String obj = array[0], catName = array[1];
					//System.out.println("Obj: \"" + obj + "\", cat: \"" + catName + "\"");
	        		if(cats.containsKey(catName)) { 
		        		Category cat = cats.get(catName);
		        		loader.insertAssignment(cat, obj);
	        		} else { 
	        			System.err.println("No category: " + catName);
	        		}
	        	}
	        	
	        	System.out.println(args[i] + " ==> " + arrays.size() + " entries.");
	        }
	        
	        loader.commitTransaction();
	        loader.close();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

    }
    
    private static Vector<String[]> loadPairs(File f) throws IOException { 
    	Vector<String[]> arrays = new Vector<String[]>();
    	BufferedReader br = new BufferedReader(new FileReader(f));
    	String line;
    	while((line = br.readLine()) != null) { 
    		line = line.trim();
    		if(line.length() > 0) { 
    			String[] array = line.split("\\s+");
				if(array.length > 2) { 
					arrays.add(array);
				}
    		}
    	}
    	br.close();
    	
    	return arrays;
    }
}
