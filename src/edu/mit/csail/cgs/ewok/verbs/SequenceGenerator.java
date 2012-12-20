package edu.mit.csail.cgs.ewok.verbs;

import java.io.IOException;
import java.io.File;
import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;

/** 
 * <code>SequenceGenerator</code> maps a Region to the genomic
 * sequence included in that Region.
 */
public class SequenceGenerator<X extends Region> implements Mapper<X,String>, SelfDescribingVerb {

    private static Map<Integer,String> cache;
    private boolean useCache = false;
    private boolean useLocalFiles = true;
    private String genomePath = null;

    private static Map<String, String[]> regionCache;
    private static Map<String, int[]> regionStarts;
    private static boolean regionIsCached = false;
    public boolean isRegionCached(){return regionIsCached;}
    
    // no longer used, but kept for compatibility 
    public SequenceGenerator (Genome g) {        
    }
    public SequenceGenerator() {}
    public void useCache(boolean b) {
        if (b && cache == null) {
            cache = new HashMap<Integer,String>();
        } 
        useCache = b;
    }
    public void useLocalFiles(boolean b) {
        useLocalFiles = b;
    }
    public void setGenomePath(String genomePath){
    	this.genomePath = genomePath;
    }
    
    /** cache the whole chromosome of this region */
    private void cache(X region) throws SQLException, IOException {
    	String chr = region.getChrom();
        int chromid = region.getGenome().getChromID(chr);
        synchronized(cache) {
            if (cache.containsKey(chromid)) {
                return;
            }
        }
        String chromseq = null;
        if (useLocalFiles) {
        	if (genomePath==null)
        		genomePath = "/scratch/" + region.getGenome().getVersion();
            File f = new File( genomePath + "/chr" + chr + ".fa");
            if (!f.exists()) {
                f = new File( genomePath+ "/chr" + chr + ".fasta");
            }
            if (f.exists()) {
                FASTAStream stream = new FASTAStream(f);
                while (stream.hasNext()) {
                    Pair<String,String> pair = stream.next();
                    String pairchrom = pair.car().replaceFirst("^chr","");
                    if (pairchrom.equals(chr)) {
                        chromseq = pair.cdr();
                        break;
                    }                    
                }
                stream.close(); 
                if (chromseq == null) {
                	System.err.println("\nchr"+chr+".fa file is not in correct FASTA format.\n");
                	System.exit(-1);
                }
            }
            else{
            	System.err.println("\nchr"+chr+".fa file is not found at "+genomePath+".\n");
            	System.exit(-1);
            }
        }
        if (chromseq == null) {
            java.sql.Connection cxn = DatabaseFactory.getConnection("core");
            PreparedStatement ps = cxn.prepareStatement("select sequence from chromsequence where id = ?");
            ps.setInt(1,chromid);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                chromseq = rs.getString(1);
            }   
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
        }
        if (chromseq == null) {
            return;
        }
        synchronized(cache) {
            if (!cache.containsKey(chromid)) {
                cache.put(chromid, chromseq);
            }
        }
    }
    /**
     * get sequence of specified region (including start and end)
     */
    public String execute(X region) {
    	if (regionIsCached)
    		return getRegionCacheSequence(region);  
    	
    	String result = null;
        String chromname = region.getChrom();
        
        try {
            Genome genome = region.getGenome();
            int chromid = genome.getChromID(chromname);
            if (useCache) {
                cache(region);
                String chromString = null;
                synchronized(cache) {
                    if (!cache.containsKey(chromid)) {
                        return null;
                    }
                    chromString = cache.get(chromid);
                }
                result = chromString.substring(region.getStart(), region.getEnd() + 1);
            }
            if (result == null) {
                java.sql.Connection cxn =
                DatabaseFactory.getConnection("core");
                PreparedStatement ps;
                int start = Math.max(region.getStart() + 1,0);
                ps = cxn.prepareStatement("select substr(sequence,?,?) from chromsequence where id = ?");
                ps.setInt(1,start);
                ps.setInt(2,region.getEnd() - region.getStart() + 1);
                ps.setInt(3,chromid);                   
                ResultSet rs = ps.executeQuery();
                if (rs.next()) {
                    result = rs.getString(1);
                } 
                rs.close();
                ps.close();
                cxn.commit();
                DatabaseFactory.freeConnection(cxn);
            }
        } catch (SQLException ex) {
            ex.printStackTrace();           
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect to core",ex);
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException("Couldn't load file to cache " + ex.toString(), ex);
        }

        if (result == null) {
            throw new DatabaseException("Couldn't get any sequence for " + region);
        }

        if (result.length() != region.getWidth()) {
            System.err.println("Wanted " + region + "(" + 
            		region.getWidth() + ") but only got " + result.length());
        }

        return result;
    }
    
    /**
     * Setup light-weight region cache of genome sequences, cover only the specified regions<br>
     * So that it does not cache the whole chromosome, save memory space. <br>
     * At the same time, retrieve some one-time sequences in rs.
     * @param regions sorted, non-overlapping regions for cache
     * @param rs regions for one-time sequence retrieval
     */
    public String[] setupRegionCache(List<Region> regions, List<Region> rs){
    	ArrayList<String> seqs = new ArrayList<String>();
    	if (regions==null||regions.isEmpty())
    		return null;    	
    	Collections.sort(regions);
    	
		// group one-time regions by chrom
		Map<String, ArrayList<Region>> chr2rs = new HashMap<String, ArrayList<Region>>();
		for(Region r:rs) {
			String chrom = r.getChrom();
			if(!chr2rs.containsKey(chrom))
				chr2rs.put(chrom, new ArrayList<Region>());
			chr2rs.get(chrom).add(r);
		}
    	
    	useCache(true);
    	regionCache = new HashMap<String, String[]>();
    	regionStarts = new HashMap<String, int[]>();
    	Genome g = regions.get(0).getGenome();
    	Region lastRegion = regions.get(regions.size()-1);
    	// setup the region space
    	String chrom = regions.get(0).getChrom();
    	int count = 0;
    	for(Region r: regions){
    		if (!r.getChrom().equals(chrom)){		// new chrom
    			regionCache.put(chrom, new String[count]);
    			regionStarts.put(chrom, new int[count]);
    			chrom = r.getChrom();
    			count = 1;
    		}
    		else		// same chrom
        		count ++;
    	}
		regionCache.put(chrom, new String[count]);
		regionStarts.put(chrom, new int[count]);
    	
		chrom = regions.get(0).getChrom();
    	count = 0;
    	for (int i=0;i<regions.size();i++){
    		Region r = regions.get(i);
    		if (!r.getChrom().equals(chrom)){	// new Chrom
    			if (cache!=null){							// for previous chrom
    				if (chr2rs.containsKey(chrom)){			// piggy-back to retrieve one-time sequences
    					for (Region r1:chr2rs.get(chrom))
    						seqs.add(execute((X)r1));
    				}
	    			synchronized(cache) {
	    				cache.put(g.getChromID(chrom), null);
	    				cache.remove(g.getChromID(chrom));	// clean cache
	    			}
	    	    	System.gc();
    			}
    			chrom = r.getChrom();
    			count = 0;
    		}		
    		// cache region using the current cache
    		synchronized(regionCache) {
    			regionStarts.get(chrom)[count]=r.getStart(); 
        		regionCache.get(chrom)[count]=execute((X)r);    			
    		}
    		count ++;
    	}
    	if (cache!=null){
			if (chr2rs.containsKey(lastRegion.getChrom())){			// piggy-back to retrieve one-time sequences
				for (Region r1:chr2rs.get(lastRegion.getChrom()))
					seqs.add(execute((X)r1));
			}
	    	synchronized(cache) {
	    		cache.put(g.getChromID(lastRegion.getChrom()), null);
	    		cache.remove(g.getChromID(lastRegion.getChrom()));
	    	}
	    	// retrieve those regions that are not in the chromosomes of cache regions
	    	for(String chr:chr2rs.keySet()){
	    		if (!regionCache.containsKey(chr)){
	    			for (Region r1:chr2rs.get(chr))
						seqs.add(execute((X)r1));
	    			synchronized(cache) {
	    				cache.put(g.getChromID(chrom), null);
	    				cache.remove(g.getChromID(chrom));	// clean cache for this chrom
	    			}
	    			System.gc();
	    		}
	    	}
	    	cache=null;
	    	System.gc();
    	}
    	
    	regionIsCached = true;
    	String[] result = new String[seqs.size()];
    	for (int i=0;i<result.length;i++)
    		result[i]=seqs.get(i);
    	return result;
    }
    
    private String getRegionCacheSequence(Region r){
    	int[] starts = regionStarts.get(r.getChrom());
    	int idx = Arrays.binarySearch(starts, r.getStart());
    	if( idx < 0 ) { idx = -idx - 2; }
    	if (!regionCache.containsKey(r.getChrom()))
    		return null;
	    synchronized(regionCache) {
	    	try{
	    		return regionCache.get(r.getChrom())[idx].substring(r.getStart()-starts[idx], r.getEnd()-starts[idx]+1);
	    	}
	    	catch(Exception e){
	    		e.printStackTrace(System.out);
	    		System.out.println(r.toString()+" idx="+idx+" starts[idx]="+starts[idx]);
	    	    return null;
	    	}
    	}
    }
    
    public static void clearCache() {
        synchronized(cache) {
            cache.clear();
        }
    }
    private static final String[] inputNames = { "Regions" };
    private static final EchoType[] inputTypes = { new ClassType(Region.class) };
    private static final EchoType outputType = new ClassType(String.class);

    public EchoType[] getInputClasses() { return inputTypes; }

    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }    
}
