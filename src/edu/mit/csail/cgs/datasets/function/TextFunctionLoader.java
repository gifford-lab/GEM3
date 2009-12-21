/*
 * Created on Nov 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.function.parsing.*;

public class TextFunctionLoader implements FunctionLoader {
    
    public static void main(String[] args) { 
        try {
            File obo = new File(args[0]);
            TextFunctionLoader loader = new TextFunctionLoader(obo);
            FunctionVersion version = loader.getVersion(obo.getName());
            
            loader.addGOAFile(new File(args[1]));
            
            System.out.print(">"); System.out.flush();
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            String line;
            while((line = br.readLine()) != null) { 
                line = line.trim();
				Category c = loader.getCategory(line);
                Collection<Assignment> assigns = loader.getAssignments(c);
                for(Assignment a : assigns) { 
                    System.out.println(String.format("\t%s --> %s (%s)",
                            a.getObject(), a.getCategory().getName(),
                            a.getCategory().getDescription()));
                }
				System.out.println(String.format("# Assigns: %d", assigns.size()));
                System.out.print(">"); System.out.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
    }
    
    private FunctionVersion version;
    private Set<Category> categories;
    private Set<Assignment> assignments;
    private Map<String,Category> categoryMap;
    private Map<Category,Set<Category>> categoryChildren;
    private Map<Category,Set<Assignment>> categoryAssignments;
    private Map<String,Set<Assignment>> objectAssignments;
    
    public TextFunctionLoader(File oboFile) throws IOException {  
        version = new FunctionVersion(oboFile.getName());
        categories = new LinkedHashSet<Category>();
        assignments = new LinkedHashSet<Assignment>();
        categoryChildren = new HashMap<Category,Set<Category>>();
        categoryAssignments = new HashMap<Category,Set<Assignment>>();
        objectAssignments = new HashMap<String,Set<Assignment>>();
        categoryMap = new HashMap<String,Category>();
        
        loadOBOFile(oboFile);
    }

	public Category getCategory(String cname) { return categoryMap.get(cname); }
	public Category getCategory(FunctionVersion fv, String cname) { return categoryMap.get(cname); }
    
    private void loadOBOFile(File f) throws IOException { 
        OBOParser parser = new OBOParser(f);
        ProvisionalGOTerm[] terms = parser.createGOTerms();
        
        for(int i = 0; i < terms.length; i++) { 
            ProvisionalGOTerm term = terms[i];
            Category c = new Category(version, term.getID(), term.getName());
            categories.add(c);
            categoryMap.put(term.getID(), c);
        }

        for(int i = 0; i < terms.length; i++) {
            ProvisionalGOTerm term = terms[i];
            Category c = categoryMap.get(term.getID());
            for(ProvisionalGOTerm pterm : term.getParents()) { 
                Category pc = categoryMap.get(pterm.getID());
                c.addParent(pc);
                
                if(!categoryChildren.containsKey(pc)) { 
                    categoryChildren.put(pc, new HashSet<Category>());
                }
                
                categoryChildren.get(pc).add(c);
            }
        }
    }
    
    public void addGOAFile(File goa) throws IOException { 
        GOAIterator itr = new GOAIterator(goa);
        while(itr.hasNext()) { 
            GOALine line = itr.next();
            Category c = categoryMap.get(line.getGOID());
            if(c != null) { 
                addAssignment(new Assignment(line.getObjectID(),c));
                addAssignment(new Assignment(line.getObjectSymbol(),c));
            } else { 
                System.err.println(String.format("Couldn't find category \"%s\"",
                        line.getGOID()));
            }
        }
    }

    public void addAssignment(Assignment a) { 
        Category c = a.getCategory();
        if(!categoryAssignments.containsKey(c)) { 
            categoryAssignments.put(c, new HashSet<Assignment>());
        }
        categoryAssignments.get(c).add(a);
        
        String obj = a.getObject();
        if(!objectAssignments.containsKey(obj)) { 
            objectAssignments.put(obj, new HashSet<Assignment>());
        }
        objectAssignments.get(obj).add(a);
        
        assignments.add(a);
    }

    public void close() {
        // do nothing.
    }
    
    public Collection<Assignment> getAllAssignments(Category c) { 
        HashSet<Assignment> assigns = new HashSet<Assignment>();
        HashSet<Category> seen = new HashSet<Category>();
        collectAllAssignments(c, assigns, seen);
        return assigns;
    }
    
    /*
     * A helper method, for getAllAssignments().
     */
    private void collectAllAssignments(Category c, 
            HashSet<Assignment> assigns, 
            HashSet<Category> catsSeen) {
        
        if(!catsSeen.contains(c)) { 
            catsSeen.add(c);
            assigns.addAll(getAssignments(c));
            Collection<Category> children = getChildCategories(c);
            for(Category child : children) { 
                collectAllAssignments(child, assigns, catsSeen);
            }
        }
    }
    
    public Collection<Category> getChildCategories(Category c) { 
        return categoryChildren.containsKey(c) ? 
                categoryChildren.get(c) : new LinkedList<Category>();
    }

    public Collection<Assignment> getAssignments(Category c) {
        return categoryAssignments.containsKey(c) ? 
            categoryAssignments.get(c) : 
            new LinkedList<Assignment>();
    }


    public Collection<FunctionVersion> getAllVersions() {
        LinkedList<FunctionVersion> vs = new LinkedList<FunctionVersion>();
        vs.add(version);
        return vs;
    }

    public Collection<Assignment> getAssignments(FunctionVersion version) {
        return assignments;
    }

    public Collection<Assignment> getAssignments(String obj, FunctionVersion fv) {
        if(fv.equals(version)) { 
            return objectAssignments.containsKey(obj) ? 
                    objectAssignments.get(obj) : 
                    new HashSet<Assignment>();
        } else { 
            return new LinkedList<Assignment>();
        }
    }

    public Collection<Category> getCategories(FunctionVersion fv) {
        return fv.equals(version) ? 
                categories : new HashSet<Category>();
    }

    public FunctionVersion getVersion(String versionName) {
        return versionName.equals(version.getName()) ? 
                version : null;
    }

}
