package edu.mit.csail.cgs.datasets.species;

// import org.biojava.ontology.*;
// import org.biojava.bio.*;
// import org.biojava.bio.symbol.*;
// import org.biojava.bio.seq.io.*;
// import org.biojava.bio.seq.*;
// import org.biojava.bio.program.gff3.*;
// import org.biojava.utils.*;
// import org.biojava.bio.seq.impl.*;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDataset;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import java.util.*;
import java.io.*;
import java.sql.*;

/**
 * A <code>Genome</code> represents one version (or genome build) of some species.
 */
public class Genome implements edu.mit.csail.cgs.utils.Closeable {
    
    public static class ChromosomeInfo { 
        
        private int length;
        private String name;
        private int dbid;
        
        private ChromosomeInfo(int id, int len, String n) { 
            length = len;
            dbid = id;
            name = n;
        }
        
        public String getName() { return name; }
        public int getDBID() { return dbid; }
        public int getLength() { return length; }
        
        public int hashCode() { 
            int code = 17;
            code += dbid; code *= 37;
            code += name.hashCode(); code *= 37;
            code += length; code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof ChromosomeInfo)) { return false; }
            ChromosomeInfo i = (ChromosomeInfo)o;
            if(dbid != i.dbid) { return false; }
            if(!name.equals(i.name)) { return false; }
            if(length != i.length) { return false; }
            return true;
        }
        
        public String toString() { return name  + " (" + length + " bp)"; }
    }
	
    private static int[] romvals;
    private static String[] intvals;
    private boolean loadedgenes = false;

    private String species, version;
    private int speciesid, dbid;
    private java.sql.Connection cxn;
    private boolean isyeast = false;

    private Map<String,ChromosomeInfo> chroms;
    private Map<Integer,ChromosomeInfo> revchroms;
    
    private ChipChipDataset dataset;
    
    /**
     * This constructor is *only* for creating a 'temporary' genome, with the 
     * given name, that does not correspond to an element in the database.
     * @param tempName
     */
    public Genome(String tempName) { 
    	species = "FakeOrganism";
    	version = tempName;
    	speciesid = dbid = -1;
    	cxn = null;
    	isyeast = false;
    	chroms = new HashMap<String,ChromosomeInfo>();
    	revchroms = new HashMap<Integer,ChromosomeInfo>();
    	
    	ChromosomeInfo info = new ChromosomeInfo(-1, 10000, "chrom");
    	chroms.put(info.getName(), info);
    	revchroms.put(-1, info);
    	
    	dataset = null;
    }

    public Genome(String tempName, Integer... lengths) { 
    	species = "FakeOrganism";
    	version = tempName;
    	speciesid = dbid = -1;
    	cxn = null;
    	isyeast = false;
    	chroms = new HashMap<String,ChromosomeInfo>();
    	revchroms = new HashMap<Integer,ChromosomeInfo>();
    	
    	for(int i = 0; i < lengths.length; i++) { 
    		ChromosomeInfo info = new ChromosomeInfo(-(i+1), lengths[i], String.format("chr%d", i+1));
        	chroms.put(info.getName(), info);
        	revchroms.put(info.dbid, info);
    	}
    	
    	dataset = null;
    }
    
    public Genome(String tempName, File chrLengths) {
    	species = "FakeOrganism";
    	version = tempName;
    	speciesid = dbid = -1;
    	cxn = null;
    	isyeast = false;
    	chroms = new HashMap<String,ChromosomeInfo>();
    	revchroms = new HashMap<Integer,ChromosomeInfo>();
    	if(!chrLengths.isFile()){System.err.println("Invalid genome info file name");System.exit(1);}
        BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(chrLengths));
		    String line;
	        int id=0;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
	            	String chr = words[0].replaceFirst("^chromosome", "");
	            	chr = chr.replaceFirst("^chrom", "");
	            	chr = chr.replaceFirst("^chr", "");
	            	ChromosomeInfo info = new ChromosomeInfo(-(id+1), Integer.parseInt(words[1]), chr);
	            	chroms.put(info.getName(), info);
	            	revchroms.put(info.dbid, info);
	            }
	    	}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	dataset = null;
    }
    
    public Genome(String tempName, Map<String, Integer> chrLengthMap) {
    	species = "FakeOrganism";
    	version = tempName;
    	speciesid = dbid = -1;
    	cxn = null;
    	isyeast = false;
    	chroms = new HashMap<String,ChromosomeInfo>();
    	revchroms = new HashMap<Integer,ChromosomeInfo>();
    	int id=0;
    	for(String s : chrLengthMap.keySet()){
    		ChromosomeInfo info = new ChromosomeInfo(id--, chrLengthMap.get(s), s);
        	chroms.put(info.getName(), info);
        	revchroms.put(info.dbid, info);
    	}
	   dataset = null;
    }
    
    public Genome(String tempSpecies, String tempVersion, Pair<String,Integer>... lengths) { 
    	species = tempSpecies;
    	version = tempVersion;
    	speciesid = dbid = -1;
    	cxn = null;
    	isyeast = false;
    	chroms = new HashMap<String,ChromosomeInfo>();
    	revchroms = new HashMap<Integer,ChromosomeInfo>();
    	
    	for(int i = 0; i < lengths.length; i++) { 
    		Pair<String,Integer> p = lengths[i];
        	ChromosomeInfo info = new ChromosomeInfo(-(i+1), lengths[i].cdr(), lengths[i].car());
        	chroms.put(info.getName(), info);
        	revchroms.put(info.dbid, info);
    	}
    	
    	dataset = null;
    }

    /**
     * Constructs a new Genome with the specified species name and genome version.
     */
    public Genome(String species, String version) throws NotFoundException {
        this.species = species;
        this.version = version;
        dataset = null; 
        chroms = null;
        revchroms = null;
        
        if (this.species.equals("Saccharomyces cerevisiae")) {
            isyeast = true;
        }
        try {
            cxn = DatabaseFactory.getConnection("core");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select id from species where name = '" + species + "'");
            if (rs.next()) {
                this.speciesid = rs.getInt(1);
            } else {
                throw new NotFoundException("Couldn't find " + species);
            }
            rs.close();
            rs = stmt.executeQuery("select id from genome where species = " + speciesid + 
                                             " and version ='" + version + "'");
            if (rs.next()) {
                dbid = rs.getInt(1);
            } else {
                throw new NotFoundException("Couldn't find " + species);
            }
            rs.close();
            stmt.close();

            fillChroms();
        
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + species + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core");
        }
    }
    
    /**
     * Constructs a new Genome from a species database identifier and a genome version.
     */
    public Genome(int speciesid, String version) throws NotFoundException {
        this.speciesid = speciesid;
        this.version = version;
        try {
            cxn = DatabaseFactory.getConnection("core");           
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select name from species where id = " + speciesid);
            if (rs.next()) {
                this.species = rs.getString(1);
            } else {
                throw new NotFoundException("Couldn't find " + species);
            }
            rs.close();
            rs = stmt.executeQuery("select id from genome where species = " + speciesid + 
                                             " and version ='" + version + "'");
            if (rs.next()) {
                dbid = rs.getInt(1);
            } else {
                throw new NotFoundException("Couldn't find " + species + " version " + version);
            }
            rs.close();
            stmt.close();

            fillChroms();
        
        } catch (SQLException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't find " + species + ": "+ ex.toString(),ex);
        }  catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect with role core");
        }

        if (this.species.equals("Saccharomyces cerevisiae")) {
            isyeast = true;
        }
    }
    /**
     * Returns a <code>ChipChipDataset</code> for this Genome.
     */
    public ChipChipDataset getChipChipDataset() { 
        if(dataset == null) { dataset = new ChipChipDataset(this); }
        return dataset;
    }

    public boolean isClosed() { 
        return cxn != null;
    }
    
    public void close() { 
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }
    
    /* retrieves the chromosomes for this Genome from the database and fills
       the relevant data structures: chroms and revchroms */
    private void fillChroms() throws SQLException {
        //        System.out.println(String.format("fillChroms() -- %s %d", version, dbid));

        chroms = new HashMap<String,ChromosomeInfo>();
        revchroms = new HashMap<Integer,ChromosomeInfo>();
        
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select c.id, c.name, length(cs.sequence) from chromosome c, chromsequence cs " +
                "where c.id=cs.id and c.genome=" + dbid);
        while(rs.next()) { 
            int dbid = rs.getInt(1);
            String name = rs.getString(2);
            int length = rs.getInt(3);
            
            ChromosomeInfo info = new ChromosomeInfo(dbid, length, name);
            if(chroms.containsKey(name) || revchroms.containsKey(dbid)) { 
                throw new IllegalArgumentException("Duplicate name \"" + name + 
                        "\" seems to exist in genome " + version);
            }
            chroms.put(name, info);
            revchroms.put(dbid, info);
            //			System.out.println(String.format("\t%s #%d (%d)", name, dbid, length));
        }
        
        rs.close();
        s.close();
    }

    public String getName() {return species; }
    public String getVersion() {return version;}
    public String getSpecies() {return species;}    
    public String getDescription() throws SQLException { 
        Statement s = cxn.createStatement();
        String desc = null;
        ResultSet rs = s.executeQuery("select description from genome where id=" + dbid);
        if(rs.next()) { desc = rs.getString(1); }
        rs.close();
        s.close();
        return desc;
    }

    /**
     * Returns the full sequence for the specified chromosome
     */
    public String getChromosomeSequence(ChromosomeInfo info) throws SQLException { 
        StringBuilder sb = new StringBuilder();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select sequence from chromsequence where id=" + info.getDBID());
        if(rs.next()) { 
            String seq = rs.getString(1);
            sb.append(seq);
        }
        rs.close();
        s.close();
        return sb.toString();        
    }
    
    /**
     * Returns a substring of the chromosomes sequence starting at <code>start</code> and of 
     * length <code>end-start+1</code>.
     */
    public String getChromosomeSequence(ChromosomeInfo info, int start, int end) throws SQLException {
        StringBuilder sb = new StringBuilder();
        Statement s = cxn.createStatement();
        
        int windowStart = start; 
        int windowLength = (end-start+1);
        
        ResultSet rs = s.executeQuery("select substr(sequence," + windowStart + "," + windowLength + ") " +
                "from chromsequence where id=" + info.getDBID());
        
        if(rs.next()) { 
            String seq = rs.getString(1);
            sb.append(seq);
        }
        
        rs.close();
        s.close();
        return sb.toString();
    }

    /** Maps a chromosome database identifier to a name */        
    public String getChromName(int chromID) { return revchroms.get(chromID).getName(); }
    /** Returns true iff this Genome contains a chromosome with the supplied database identifier */
    public boolean containsChromID(int chromID) { return revchroms.containsKey(chromID); }
    /** Returns true iff this Genome contains a chromosome with the supplied name */
    public boolean containsChromName(String chromName) { return chroms.containsKey(chromName); }
    /** Maps a chromosome name to its database identifier */
    public int getChromID(String chromName) { 
        if (chroms.get(chromName) == null) {
            throw new NullPointerException("Null chromosome for " + chromName);
        }
        return chroms.get(chromName).getDBID(); 
    }
    /** Returns the complete mapping from chromosome names to DBIDs for this Genome */
    public Map<String,Integer> getChromIDMap() { 
        Map<String,Integer> chromID = new HashMap<String,Integer>();
        for(String n : chroms.keySet()) { chromID.put(n, chroms.get(n).getDBID()); }
        return chromID;
    }
    /** Returns the complete mapping from chromosome DBIDs to names for this Genome */
    public Map<Integer,String> getRevChromIDMap() { 
        Map<Integer,String> chromID = new HashMap<Integer,String>();
        for(int dbid : revchroms.keySet()) { chromID.put(dbid, revchroms.get(dbid).getName()); }
        return chromID;
    }
    /** Returns the length of the specified chromosome, in base pairs */
    public int getChromLength(String chromName) { return chroms.get(chromName).getLength(); }
    /** Returns the complete map from chromosome name to length */
    public Map<String,Integer> getChromLengthMap() { 
        Map<String,Integer> chromLengths = new HashMap<String,Integer>();
        for(String n : chroms.keySet()) { chromLengths.put(n, chroms.get(n).getLength()); }
        return chromLengths;
    }
    /** Returns the genome info string with chromosome name <tab> length format*/
    public String getGenomeInfo(){
    	StringBuilder sb = new StringBuilder();
    	for(String n : chroms.keySet()) { sb.append(n).append("\t").append(chroms.get(n).getLength()).append("\n"); }
        return sb.toString();
    }
    /** Converts a goofy chromosome name (eg, chrIV) to a more sensible name (eg 4) */
    public String fixChromName(String chrom) {return fixChrom(chrom);}
    /** Maps a sensible chromosome name back to the roman numeral form */
    public String unfixChromName(String chrom) {return unfixChrom(chrom);}
    public String unfixChromName(int chrom) {return unfixChrom(chrom);}
    
    public double getGenomeLength() { 
        double totalLen=0;
        for(String n : chroms.keySet()) { totalLen+= (double)chroms.get(n).getLength();}
        return totalLen;
    }


    public static String unfixChrom(String c) {
        return unfixChrom(Integer.parseInt(c));
    }

    public static String unfixChrom(int chrom) {
        if(intvals == null) { 
            intvals = new String[10];
            intvals[0] = "X";
            intvals[1] = "I";
            intvals[2] = "II";
            intvals[3] = "III";
            intvals[4] = "IV";
            intvals[5] = "V";
            intvals[6] = "VI";
            intvals[7] = "VII";
            intvals[8] = "VIII";
            intvals[9] = "IX";
        }
        
        StringBuilder sb = new StringBuilder();
        sb.append("chr");
        
        while(chrom >= 10) { 
            chrom -= 10;
            sb.append(intvals[0]);
        }
        
        if(chrom > 0) { 
            sb.append(intvals[chrom]);
        }
        
        return sb.toString();
    }
    
    public static String fixChrom(String chrom) {
        if (romvals == null) {
            romvals = new int[Character.getNumericValue('Z')];
            romvals[Character.getNumericValue('X')] = 10;
            romvals[Character.getNumericValue('V')] = 5;
            romvals[Character.getNumericValue('I')] = 1;
        }
        String chr = chrom;
        chr = chr.replaceAll("\\.fa?s?t?a$","");
        if (chr.matches("^[cC][hH][rR].*")) {
            chr = chr.substring(3);
        } 
//         if (isyeast && chr.matches("^[XVI]+$")) {
//             int val = 0, pos = 1, curval, lastval, buffer; char cur, last;
//             boolean random = false;
//             if (chr.matches("_random$")) {
//                 random = true;
//                 chr.replaceFirst("_random$","");
//             }            
//             last = chr.charAt(0);
//             lastval = romvals[Character.getNumericValue(last)];
//             buffer = lastval;
//             //            System.err.println("== " + buffer);
//             while (pos < chr.length()) {
//                 cur = chr.charAt(pos);
//                 curval = romvals[Character.getNumericValue(cur)];
//                 if (curval > lastval) {
//                     val += curval - lastval;
//                     buffer = 0;
//                 } else if (cur != last) {
//                     val += buffer;
//                     buffer = curval;
//                 } else {
//                     buffer += curval;
//                 }
//                 //                System.err.println(pos + ":" + cur + "," + curval + "," + buffer + "," + val);
//                 last = cur;
//                 lastval = curval;
//                 pos++;
//             }
//             val += buffer;
//             //            System.err.println(",,, " + buffer + "," + val);
//             if (random) {
//                 return Integer.toString(val) + "_random";
//             } else {
//                 return Integer.toString(val);
//             }
//         } else 
        if (chr.matches("^[1234567890MmtUnXY]+(_random)?[LRh]?$")) {
            return chr;
        } else {
            throw new NumberFormatException("Can't fix chrom name " + chrom + "," + chr);
        }
    }

    public static String fixYeastChrom(String chrom) {
        if (romvals == null) {
            romvals = new int[Character.getNumericValue('Z')];
            romvals[Character.getNumericValue('X')] = 10;
            romvals[Character.getNumericValue('V')] = 5;
            romvals[Character.getNumericValue('I')] = 1;
        }
        String chr = chrom;
        chr = chr.replaceAll("\\.fa?s?t?a$","");
        if (chr.matches("^[cC][hH][rR].*")) {
            chr = chr.substring(3);
        } 
        if (chr.matches("^[XVI]+$")) {
            int val = 0, pos = 1, curval, lastval, buffer; char cur, last;
            boolean random = false;
            if (chr.matches("_random$")) {
                random = true;
                chr.replaceFirst("_random$","");
            }            
            last = chr.charAt(0);
            lastval = romvals[Character.getNumericValue(last)];
            buffer = lastval;
            //            System.err.println("== " + buffer);
            while (pos < chr.length()) {
                cur = chr.charAt(pos);
                curval = romvals[Character.getNumericValue(cur)];
                if (curval > lastval) {
                    val += curval - lastval;
                    buffer = 0;
                } else if (cur != last) {
                    val += buffer;
                    buffer = curval;
                } else {
                    buffer += curval;
                }
                //                System.err.println(pos + ":" + cur + "," + curval + "," + buffer + "," + val);
                last = cur;
                lastval = curval;
                pos++;
            }
            val += buffer;
            //            System.err.println(",,, " + buffer + "," + val);
            if (random) {
                return Integer.toString(val) + "_random";
            } else {
                return Integer.toString(val);
            }
        } else 
        if (chr.matches("^[1234567890MUXY]+(_random)?[LRh]?$")) {
            return chr;
        } else if (chr.matches("Mito")) {
            return "mt";
        } else {
            throw new NumberFormatException("Can't fix chrom name " + chrom + "," + chr);
        }
    }

    /** Returns a list of all chromosome names in this GEnome */
    public List<String> getChromList() { return new LinkedList<String>(chroms.keySet()); }
    /** Maps a chromosome name to the corresponding <code>ChromosomeInfo</code> object */
    public ChromosomeInfo getChrom(String name) { return chroms.get(name); }
    /** Returns the database identifier for this Genome */
    public int getDBID() {return dbid;}
    public int getSpeciesDBID() { return speciesid; }
    
    /** Returns a read connection to the UCSC database for this
     * genome 
     */
    public java.sql.Connection getUcscConnection() throws SQLException {
        try {
        	if(getVersion().equals("vv1")) { 
        		return DatabaseFactory.getConnection("vaccinia_vv1");
        	}
        	
            String v = this.getVersion().replaceAll("[^\\w\\-]","_");
            return DatabaseFactory.getConnection("ucsc_" + v);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't create a database connection for genome " + 
                                        getVersion(),ex);
        }

    }
    
    public String toString() {
        return getSpecies() +","+getVersion();
    }
    
    public int hashCode() {
        return getSpecies().hashCode()*37 + getVersion().hashCode();
    }

    public boolean equals(Object o) {
        if (o instanceof Genome) {
            Genome other = (Genome)o;
            return (getSpecies().equals(other.getSpecies()) &&
                    getVersion().equals(other.getVersion()));
        } else {
            return false;
        }
    }
}

