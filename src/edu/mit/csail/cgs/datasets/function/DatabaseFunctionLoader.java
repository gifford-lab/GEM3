/**
 * 
 */
package edu.mit.csail.cgs.datasets.function;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.io.*;

import edu.mit.csail.cgs.datasets.function.parsing.MIPSTextFunctionalCategories;
import edu.mit.csail.cgs.datasets.function.parsing.OBOParser;
import edu.mit.csail.cgs.utils.Enrichment;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author Timothy Danford
 */
public class DatabaseFunctionLoader 
    implements edu.mit.csail.cgs.utils.Closeable, FunctionLoader {
    
    public static String dbRole;
    
    static { 
        dbRole = "function";
    }
    
    public static void main(String[] args) {
        try {
            DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
            FunctionalUtils utils = new FunctionalUtils(loader, "LiverAnnotations");
            
            BufferedReader br = new BufferedReader(new FileReader(new File(args[0])));
            
            Set<String> total = new HashSet<String>();
            Set<String> target = new HashSet<String>();
            
            String line = null;
            while((line = br.readLine()) != null) { 
            	line = line.trim();
            	if(line.length() > 0) { 
            		String[] array = line.split("\\s+");
            		boolean included = array[0].equals("1");
            		String val = array[1];
            		total.add(val);
            		if(included) { target.add(val); }
            	}
            }
            
            br.close();
            
            Map<String,Enrichment> enrs = utils.calculateTotalEnrichments(total, target);
            TreeSet<Enrichment> sorted = new TreeSet<Enrichment>(enrs.values());
            double logThresh = Math.log(0.01);
            
            for(Enrichment e : sorted) { 
            	if(e.getLogPValue() <= logThresh) { 
            		System.out.println(e);
            	}
            }
        
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (IOException e) {
			e.printStackTrace();
		}
        
    }

    /**
     * func_assignment:
     * id number(10)
     * version number(10)
     * object varchar(100)
     * category number(10)
     */
	
    /**
     * func_category:
     * id number(10)
     * version number(10)
     * name varchar (100)
     * description varchar(1000)
     */
	
    /**
     * func_subcategory:
     * child_id number(10)
     * parent_id number(10)
     * version number(10)
     */
	
    /**
     * func_version:
     * id number(10)
     * name varchar(500)
     */

    private Map<Integer,FunctionVersion> id2Version;
    private Map<Integer,Collection<Category>> versionID2Categories;
    private Map<Integer,Category> id2Category;
    private Map<Integer,Collection<Assignment>> categoryID2Assignments;
	
    private java.sql.Connection cxn;
	
    private PreparedStatement loadAllVersions, loadVersionByID;
    private PreparedStatement loadCategByID, loadCategsByVersion, loadCategsByNameVersion, loadCategParentsByID, loadCategChildrenByID;
    private PreparedStatement loadAssignsByCategory, loadAssignsByObjVersion, loadAssignsByVersion;
    
    private PreparedStatement insertAssign, insertCategory, insertCategoryParent, insertVersion;
	
    public DatabaseFunctionLoader() throws SQLException, UnknownRoleException {
        cxn = DatabaseFactory.getConnection(dbRole);
        init();
    }
    
    private void init() throws SQLException {   
        id2Version = new HashMap<Integer,FunctionVersion>();
        versionID2Categories = new HashMap<Integer,Collection<Category>>();
        id2Category = new HashMap<Integer,Category>();
        categoryID2Assignments = new HashMap<Integer,Collection<Assignment>>();
		
        loadAllVersions = FunctionVersion.prepareReadAll(cxn);
        loadVersionByID = FunctionVersion.prepareReadByID(cxn);
        loadCategByID = Category.prepareReadByID(cxn);
        loadCategsByNameVersion = Category.prepareReadByNameVersion(cxn);
        loadCategsByVersion = Category.prepareReadIDByVersion(cxn);
        loadCategParentsByID = Category.prepareReadParentsByID(cxn);
        loadCategChildrenByID = Category.prepareReadChildrenByID(cxn);
        loadAssignsByCategory = Assignment.createReadByCategory(cxn);
        loadAssignsByObjVersion = Assignment.createReadIDByObjectVersion(cxn);
        loadAssignsByVersion = Assignment.createReadByVersion(cxn);
        
        insertAssign = insertCategory = insertCategoryParent = insertVersion = null;
        insertAssign = Assignment.prepareInsert(cxn);
        insertCategory = Category.prepareInsert(cxn);
        insertCategoryParent = Category.prepareInsertParents(cxn);
        insertVersion = FunctionVersion.prepareInsertStatement(cxn);
	}
    
    public void close() { 
        try {
			
            loadAllVersions.close(); loadAllVersions = null;
            loadVersionByID.close(); loadVersionByID = null;
            loadCategByID.close(); loadCategByID = null;
            loadCategsByVersion.close(); loadCategsByVersion = null;
            loadCategParentsByID.close(); loadCategParentsByID = null;
            loadCategChildrenByID.close(); loadCategChildrenByID = null;
            loadCategsByNameVersion.close(); loadCategsByNameVersion = null;
            loadAssignsByCategory.close(); loadAssignsByCategory = null;
            loadAssignsByObjVersion.close(); loadAssignsByObjVersion = null;
            loadAssignsByVersion.close(); loadAssignsByVersion = null;
            
			insertAssign.close(); insertAssign = null;
			insertCategory.close(); insertCategory = null;
			insertCategoryParent.close(); insertCategoryParent = null;
			insertVersion.close(); insertVersion = null;

			DatabaseFactory.freeConnection(cxn);
			
        } catch (SQLException e) {
            e.printStackTrace();
        }
		
        cxn = null;
    }
	
    public boolean isClosed() { 
        return cxn == null;
    }
    
    /**
     * Deletes a given FunctionVersion (and all its associated data) from the database.
     * @param fv  The FunctionVersion to delete.
     * @throws SQLException
     */
    public void deleteFunctionVersion(FunctionVersion fv) throws SQLException {
        if(fv.getID() == -1) { throw new IllegalArgumentException(); }
        System.out.println("Deleting FunctionVersion \"" + fv.getName() + "\"");
        boolean oldCommit = cxn.getAutoCommit();
        cxn.setAutoCommit(false);        
        Statement s = cxn.createStatement();
        
        System.out.println("\tDeleting Assignments...");
        s.executeUpdate("delete from func_assignment where version=" + fv.getID());
        
        System.out.println("\tDeleting Category Relationships...");
        s.executeUpdate("delete from func_subcategory where version=" + fv.getID());
        
        System.out.println("\tDeleting Categories...");
        s.executeUpdate("delete from func_category where version=" + fv.getID());
        
        System.out.println("\tDeleting FunctionVersion entry...");
        s.executeUpdate("delete from func_version where id="+ fv.getID());
        
        s.close();
        cxn.commit();
        cxn.setAutoCommit(oldCommit);        
        System.out.println("Done.");
    }
	
    /**
     * Returns the list of *all* FunctionVersions that are available in the database.
     * 
     * @return A collection of FunctionVersion objects.
     * @see edu.mit.csail.cgs.datasets.function.FunctionLoader#getAllVersions()
     */
    public Collection<FunctionVersion> getAllVersions() throws SQLException { 
        if(id2Version.size() == 0) {
            try { 
                synchronized (loadAllVersions) {
                    ResultSet rs = loadAllVersions.executeQuery();
                    while(rs.next()) { 
                        FunctionVersion fv = new FunctionVersion(rs);
                        if(fv != null) { getCategories(fv); } 
                        id2Version.put(fv.getID(), fv);
                    }
                    rs.close();
                }
            } catch(Exception se) { 
                id2Version.clear();
                //throw se;
				se.printStackTrace(System.err);
            }
        }
        return id2Version.values();
    }

    /**
     * Load a FunctionVersion object by database id.
     * @param vid The DBID of the function version.
     * @return The FunctionVersion object.
     * @throws SQLException
     */
    public FunctionVersion getVersion(int vid) throws SQLException {
		FunctionVersion fv = null;
        synchronized (loadVersionByID) {
            loadVersionByID.setInt(1, vid);
            ResultSet rs = loadVersionByID.executeQuery();
            if(rs.next()) { 
                fv = new FunctionVersion(rs);
                id2Version.put(fv.getID(), fv);
            }
            if(fv != null) { getCategories(fv); } 
            rs.close();
        }
        return fv;
    }
	
    /**
     * Load a FunctionVersion object by name.
     * @param n The name of the function version.
     * @return The FunctionVersion object.
     * @throws SQLException
     */
    public FunctionVersion getVersion(String n) throws SQLException { 
        PreparedStatement ps = cxn.prepareStatement("select id from func_version where name = ?");
        ps.setString(1,n);
    	ResultSet rs = ps.executeQuery();
    	FunctionVersion fv = null;
    	if(rs.next()) { 
    		int id = rs.getInt(1);
    		fv = getVersion(id);
    	}
    	
    	rs.close();
    	ps.close();
    	
    	if(fv != null) { getCategories(fv); } 
    	
    	return fv;
    }
    
    public void insertVersion(String n) throws SQLException {
        synchronized (insertVersion) {
            insertVersion.setString(1, n);
            insertVersion.executeUpdate();
        }
    }
    
    /**
     * Returns all the Category objects associated with a particular FunctionVersion.
     * 
     * @param fv The name of the version.
     * @see edu.mit.csail.cgs.datasets.function.FunctionLoader#getCategories(edu.mit.csail.cgs.datasets.function.FunctionVersion)
     */
    public Collection<Category> getCategories(FunctionVersion fv) throws SQLException {
        if(!versionID2Categories.containsKey(fv.getID())) {
            LinkedList<Category> lst = new LinkedList<Category>();
            versionID2Categories.put(fv.getID(), lst);

            synchronized (loadCategsByVersion) {
                loadCategsByVersion.setInt(1, fv.getID());
                ResultSet rs = loadCategsByVersion.executeQuery();
                System.out.println("Loading categories (#" + fv.getID() + ")...");

                while(rs.next()) { 
                    Category c = new Category(rs, this);
                    id2Category.put(c.getID(), c);
                    lst.addLast(c);
                }
                rs.close();
            }
            
            System.out.println("Loading category relationships...");
            Statement s = cxn.createStatement();
            ResultSet prs = s.executeQuery("select child_id, parent_id from func_subcategory");
            while(prs.next()) { 
                int ci = prs.getInt(1), pi = prs.getInt(2);
                Category cp = id2Category.get(pi);
                cp.addParent(id2Category.get(ci));
            }
            
            prs.close();
            s.close();
			System.out.println("Done.");
        }

        return versionID2Categories.get(fv.getID());
    }
    
    Category getCategory(int id) { 
    	return id2Category.get(id);
    }
	
    public Category getCategory(FunctionVersion fv, String n) throws SQLException { 
    	for(Category c : versionID2Categories.get(fv.getID())) { 
    		if(c.getName().equals(n)) { 
    			return c;
    		}
    	}
    	return null;
    }
    
    public void insertCategoryRelationship(Category child, Category parent) throws SQLException { 
        synchronized (insertCategoryParent) {
            insertCategoryParent.setInt(1, child.getID());
            insertCategoryParent.setInt(2, parent.getID());
            insertCategoryParent.executeUpdate();
        }
    }

    public void insertCategory(FunctionVersion version, String name, String description) throws SQLException { 
        synchronized (insertCategory) {
            insertCategory.setInt(1, version.getID());
            insertCategory.setString(2, name);
            insertCategory.setString(3, description);
            insertCategory.executeUpdate();
        }
    }
   
    /**
     * Collects the set of *all* assignment objects which are assigned to either the 
     * given category, or *any* descendant of this category.
     * 
     * @param c  The given category.
     * @return The set of Assignment objects.
     * @throws SQLException
     */
    public Collection<Assignment> getAllAssignments(Category c) throws SQLException { 
        HashSet<Assignment> assigns = new HashSet<Assignment>();
        HashSet<Category> seen = new HashSet<Category>();
        collectAllAssignments(c, assigns, seen);
        return assigns;
    }
    
    public void insertAssignment(Category c, String obj) throws SQLException { 
        synchronized (insertAssign) {
            insertAssign.setString(1, obj);
            insertAssign.setInt(2, c.getID());
            insertAssign.executeUpdate();
        }
    }
    
    public void beginTransaction() throws SQLException { 
    	cxn.setAutoCommit(false);
    }
    
    public void commitTransaction() throws SQLException { 
    	cxn.commit();
    	cxn.setAutoCommit(true);
    }
    
    /**
     * Loads all the Assignment objects associated with this particular version.
     * 
     * @param fv  The indicated FunctionVersion
     * @return A collection of Assignment objects loaded.
     * @throws SQLException
     */
    public Collection<Assignment> getAssignments(FunctionVersion fv) throws SQLException { 
        HashSet<Assignment> assigns = new HashSet<Assignment>();
        synchronized (loadAssignsByVersion) {
            loadAssignsByVersion.setInt(1, fv.getID());
            ResultSet rs = loadAssignsByVersion.executeQuery();
            while(rs.next()) {   
                Assignment assign = new Assignment(rs, this);
                assigns.add(assign);
            }
            rs.close();
        }
        return assigns;
    }

    /*
     * A helper method, for getAllAssignments().
     */
    private void collectAllAssignments(Category c, 
                                       HashSet<Assignment> assigns, HashSet<Category> catsSeen) throws SQLException {
        
        if(!catsSeen.contains(c)) { 
            catsSeen.add(c);
            assigns.addAll(getAssignments(c));
            Collection<Category> children = getChildCategories(c);
            for(Category child : children) { 
                collectAllAssignments(child, assigns, catsSeen);
            }
        }
    }
    
    /**
     * Returns the Category objects which are the (immediate) children of the given Category.
     * @param c The given category.
     * @return A set of Category objects.
     * @throws SQLException
     */
    public Collection<Category> getChildCategories(Category c) throws SQLException { 
        HashSet<Category> cats = new HashSet<Category>();
        synchronized (loadCategChildrenByID) {
            loadCategChildrenByID.setInt(1, c.getID());
            ResultSet rs = loadCategChildrenByID.executeQuery();
            while(rs.next()) { 
                int cid = rs.getInt(1);
                Category cat = getCategory(cid);
                cats.add(cat);
            }
            rs.close();
        }
        return cats;
    }
	
    /**
     * Returns all the Assignment objects which are *immediately* associated with the 
     * given Category. 
     * 
     * @param c The given Category.
     * @return The set of Assignment objects.
     * @throws SQLException
     */
    public Collection<Assignment> getAssignments(Category c) 
        throws SQLException {
        
        if(c.getID() == -1) { throw new IllegalArgumentException(); }
        if(!categoryID2Assignments.containsKey(c.getID())) { 
            LinkedList<Assignment> assigns = new LinkedList<Assignment>();
			synchronized (loadAssignsByCategory) {
                loadAssignsByCategory.setInt(1, c.getID());
                ResultSet rs = loadAssignsByCategory.executeQuery();
                try { 
                while(rs.next()) { 
                    Assignment a = new Assignment(rs, this);
                    assigns.addLast(a);
                }
                } catch(SQLException se) { 
                    throw se;
                } finally { 
                    rs.close();
                }
            }
			
            categoryID2Assignments.put(c.getID(), assigns);
        }
        return categoryID2Assignments.get(c.getID());
    }

    /**
     * Returns the set of all Assignment objects which correspond to a particular object.
     * In essence, this is used to get you all the Category objects to which a particular element
     * is assigned.
     * 
     * @param obj The name of the object.
     * @param fv  The name of the version to search.
     * @return The set of Assignment objects, whose "Object" field is equal to obj
     */
    public Collection<Assignment> getAssignments(String obj, FunctionVersion fv) 
        throws SQLException {
		
        //System.out.println("Searching for assignments \"" + obj + "\"");
        if(fv.getID() == -1) { throw new IllegalArgumentException(); }
        LinkedList<Assignment> assigns = new LinkedList<Assignment>();
		
        synchronized (loadAssignsByObjVersion) {
            loadAssignsByObjVersion.setString(1, obj);
            loadAssignsByObjVersion.setInt(2, fv.getID());
            
            ResultSet rs = loadAssignsByObjVersion.executeQuery();
            try { 
                while(rs.next()) { 
                    Assignment a = new Assignment(rs, this);
                    assigns.addLast(a);
                }
            } catch(SQLException se) { 
                throw se;
            } finally { 
                rs.close();
            }
        }
		
        return assigns;
    }
    
    private static class CategoryComparator implements Comparator<Category> {
        
        public CategoryComparator() {}
        
        public int compare(Category c1, Category c2) { 
            if(c1.isChildCategory(c2)) { return 1; }
            if(c2.isChildCategory(c1)) { return -1; }
            return c1.getName().compareTo(c2.getName());
        }
    }
}
