/*
 * Created on Dec 20, 2005
 */
package edu.mit.csail.cgs.datasets.function;

import java.util.*;
import java.text.*;
import java.sql.*;
import java.io.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.probability.Hypergeometric;
import edu.mit.csail.cgs.datasets.function.*;

/**
 * @author tdanford
 */
public class FunctionalUtils implements edu.mit.csail.cgs.utils.Closeable {
    
    public static void main(String[] args) {
        File test = new File(args[0]);
        File bg = new File(args[1]);
        String version = args.length > 2 ? args[2] : "GO:5_23_2006";
        Set<String> ids = new HashSet<String>();
        Set<String> total = new HashSet<String>();
        try { 
            BufferedReader br = new BufferedReader(new FileReader(bg));
            String line;
            while((line = br.readLine()) != null) { 
                line = line.trim();
                if(line.length() > 0) { 
                    total.add(line);
                }
            }
            br.close();
            
            br = new BufferedReader(new FileReader(test));
            while((line = br.readLine()) != null) { 
                line = line.trim();
                if(line.length() > 0 && total.contains(line)) { 
                    ids.add(line);
                }
            }
            
			FunctionLoader loader = null;
			if(args.length > 3) { 
				File f= new File(args[2]);
				TextFunctionLoader tfl = new TextFunctionLoader(f);
				tfl.addGOAFile(new File(args[3]));
				loader = tfl;
				version = f.getName();
				//System.out.println("loading Text FunctionLoader");
			} else { 
				loader = new DatabaseFunctionLoader();	
				//System.out.println("loading DB FunctionLoader");
			}

            FunctionalUtils utils = new FunctionalUtils(loader, version);
            
            Map<String,Enrichment> enrichs = utils.calculateTotalEnrichments(total, ids);
            TreeSet<Enrichment> sorted = new TreeSet<Enrichment>(enrichs.values());
            
            for(Enrichment e : sorted) { 
                System.out.println(e.toString());
            }
            
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }
    }
	
    private FunctionLoader loader;
    private Map<String,FunctionVersion> versionMap;
    private Hypergeometric hg;
    
    private FunctionVersion version;
    private Map<String,Category> categories;
    private SetTools<String> stringTools;

    public FunctionalUtils(FunctionLoader l, String versionName) throws SQLException {
        loader = l;
        versionMap = new HashMap<String,FunctionVersion>();
        hg = new Hypergeometric();
        stringTools = new SetTools<String>();
        
        //for(FunctionVersion fv : loader.getAllVersions()) { 
            //versionMap.put(fv.getName(), fv);
        //}
		versionMap.put(versionName, loader.getVersion(versionName));
        System.err.println("Creating FU for " + versionName + " as " + versionMap.get(versionName));
        setVersion(versionName);
    }
    
    public void setVersion(String v) throws SQLException {
        if(!versionMap.containsKey(v)) { throw new IllegalArgumentException(v); }
        version = versionMap.get(v);
        System.err.println("VersionMap gave " + version);
        categories = new HashMap<String,Category>();
        for(Category c : loader.getCategories(version)) { 
            categories.put(c.getName(), c);
        }
        System.err.println("Added " + categories.size() + " categories to FunctionalUtils.");
    }
    
    public void createNewVersionFromCategories(String newVersionName, FunctionVersion oldVersion) throws SQLException  { 
    	Map<Integer,Integer> old2new = new HashMap<Integer,Integer>();
    }
    
    public FunctionLoader getFunctionLoader() { return loader; }
    
    public Category getCategory(String cname) { return categories.get(cname); }
    
    public static Set<String> convertAssignmentsToObjects(Collection<Assignment> assigns) { 
        HashSet<String> objs = new HashSet<String>();
        for(Assignment a : assigns) { objs.add(a.getObject()); }
        return objs;
    }
    
    public static Set<String> convertCategoriesToNames(Collection<Category> cats) { 
        HashSet<String> catNames = new HashSet<String>();
        for(Category c : cats) { catNames.add(c.getName()); }
        return catNames;        
    }
    
    public static int sumCountMap(Map<String,Integer> countMap) { 
        int c = 0;
        for(String k : countMap.keySet()) { c += countMap.get(k); }
        return c;
    }
    
    public Map<String,Set<String>> getImmediateCategories(Set<String> objs) throws SQLException {
        Map<String,Set<String>> assignments = new HashMap<String,Set<String>>();
        for(String o : objs) { 
            for(Assignment a : loader.getAssignments(o, version)) { 
                Category c = a.getCategory();
                if(!assignments.containsKey(c)) { 
                    assignments.put(c.getName(), new HashSet<String>());
                }
                assignments.get(c.getName()).add(a.getObject());
            }
        }
        return assignments;
    }
    
    public Map<String,Set<String>> getTotalCategories(Set<String> objs) throws SQLException {
        Map<String,Set<String>> assignments = new HashMap<String,Set<String>>();
        for(String o : objs) { 
            for(Assignment a : loader.getAssignments(o, version)) { 
                Category c = a.getCategory();
                recursivelyAssign(assignments, a.getObject(), c);
            }
        }
        return assignments;
    }
    
    private void recursivelyAssign(Map<String,Set<String>> assignments, String obj, Category c) { 
        String cname = c.getName();
        if(!assignments.containsKey(cname)) { assignments.put(cname, new HashSet<String>()); }
        if(!assignments.get(cname).contains(obj)) { 
            assignments.get(cname).add(obj);
            for(Category p : c.getParents()) { 
                recursivelyAssign(assignments, obj, p);
            }
        }
    }
    
    public Map<String,Integer> getCountedImmediateCategories(Set<String> objs) throws SQLException {
        Map<String,Integer> counts = new HashMap<String,Integer>();
        Map<String,Set<String>> assignments = getImmediateCategories(objs);
        for(String k : assignments.keySet()) { 
            counts.put(k, assignments.get(k).size());
        }
        return counts;
    }
    
    public Map<String,Integer> getCountedTotalCategories(Set<String> objs) throws SQLException {    	
        Map<String,Integer> counts = new HashMap<String,Integer>();
        Map<String,Set<String>> assignments = getTotalCategories(objs);
        for(String k : assignments.keySet()) { 
            counts.put(k, assignments.get(k).size());
        }
        return counts;
    }
    
    public Map<String,Enrichment> calculateImmediateEnrichments(Set<String> total, Set<String> objs) 
        throws SQLException { 
        
        Set<String> filteredObjs = stringTools.intersection(total, objs);
        
        Map<String,Integer> totalCounts = getCountedImmediateCategories(total);
        Map<String,Integer> objCounts = getCountedImmediateCategories(filteredObjs);
    	Map<String,Enrichment> pvalues = new HashMap<String,Enrichment>();
        
        int N = total.size();
        int n = filteredObjs.size();
        
        for(String cname : totalCounts.keySet()) { 
            int theta = totalCounts.get(cname);
            int x = objCounts.containsKey(cname) ? objCounts.get(cname) : 0;
            double logpv = hg.log_hypgeomPValue(N, theta, n, x);
            Category c = getCategory(cname);
    		pvalues.put(cname, new Enrichment(c.getName() + "\t" + c.getDescription(), N, theta, n, x, logpv));
        }
        
        return pvalues;
    }
    
    public Map<String,Enrichment> calculateTotalEnrichments(Set<String> total, Set<String> objs) 
        throws SQLException { 

        Set<String> filteredObjs = stringTools.intersection(total, objs);

        int N = total.size();
        int n = filteredObjs.size();
		//System.out.println(String.format("Total %d, Filtered Test %d", N, n));
        
        Map<String,Integer> totalCounts = getCountedTotalCategories(total);
        Map<String,Integer> objCounts = getCountedTotalCategories(filteredObjs);
        Map<String,Enrichment> pvalues = new HashMap<String,Enrichment>();

		/*
		for(String k : objCounts.keySet()) { 
			System.out.println(String.format("%s --> %d,%d", k, objCounts.get(k), totalCounts.get(k)));
		}
		*/
        
        for(String cname : totalCounts.keySet()) { 
            int theta = totalCounts.get(cname);
            int x = objCounts.containsKey(cname) ? objCounts.get(cname) : 0;
            double logpv = hg.log_hypgeomPValue(N, theta, n, x);
            Category c = getCategory(cname);
			if(c != null) { 
				pvalues.put(cname, new Enrichment(c.getName() + "\t" + c.getDescription(), N, theta, n, x, logpv));
			} else { 
				//System.err.println("ERROR!  Couldn't find category for name " + cname);
			}
        }
        
        return pvalues;
    }
    
    public int countImmediateOverlap(Set<String> objs, Category c) throws SQLException { 
        int overlap = 0;
        for(Assignment a : loader.getAssignments(c)) { 
            if(objs.contains(a.getObject())) { 
                overlap += 1; 
            }
        }
        return overlap;
    }

    public int countTotalOverlap(Set<String> objs, Category c) throws SQLException { 
        int overlap = 0;
        for(Assignment a : loader.getAllAssignments(c)) { 
            if(objs.contains(a.getObject())) { 
                overlap += 1; 
            }
        }
        return overlap;
    }

    public Enrichment calculateImmediateLogEnrichment(Set<String> total, Set<String> objs, Category c) 
        throws SQLException { 
        
        Set<String> filteredObjs = stringTools.intersection(total, objs);

        int N = 0;
        int n = 0;
        for(Assignment ta : loader.getAssignments(version)) { 
            if(total.contains(ta.getObject())) { N += 1; }
            if(filteredObjs.contains(ta.getObject())) { n += 1; }
        }
        
        int theta = countImmediateOverlap(total, c);
        int x = countImmediateOverlap(filteredObjs, c);

        double logpv = hg.log_hypgeomPValue(N, theta, n, x);
        return new Enrichment(c.getName() + "\t" + c.getDescription(), N, theta, n, x, logpv);
    }
    
    public Enrichment calculateTotalLogEnrichment(Set<String> total, Set<String> objs, Category c) 
        throws SQLException { 
        
        Set<String> filteredObjs = stringTools.intersection(total, objs);
        
        int N = 0;
        int n = 0;
        for(Assignment ta : loader.getAssignments(version)) { 
            if(total.contains(ta.getObject())) { N += 1; }
            if(filteredObjs.contains(ta.getObject())) { n += 1; }
        }
        
        int theta = countTotalOverlap(total, c);
        int x = countTotalOverlap(filteredObjs, c);
        
        double logpv = hg.log_hypgeomPValue(N, theta, n, x);
        return new Enrichment(c.getName() + "\t" + c.getDescription(), N, theta, n, x, logpv);
    }

    private static DecimalFormat nf;
    
    static { 
    	nf = (DecimalFormat)DecimalFormat.getInstance();
    	nf.setMaximumIntegerDigits(1);
    	nf.setMinimumIntegerDigits(1);
    	nf.setMaximumFractionDigits(6);
    	nf.setMinimumFractionDigits(6);
    }    
    public void close() {
        loader.close();
        loader = null;
    }
    public boolean isClosed() {
        return (loader == null);
    }
}
