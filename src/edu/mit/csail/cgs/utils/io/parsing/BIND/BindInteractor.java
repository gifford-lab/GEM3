/*
 * Created on Sep 7, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.BIND;

import java.util.*;

/**
 * @author tdanford
 */
public class BindInteractor {
    
    private int id;
    private String accession;
    private String db;
    private int taxon;
    private String type;
    
    private Set<String> shortNames, otherNames;

    public BindInteractor(int id, String acc, String db, String type, int taxon) {
        this.id = id;
        accession = acc;
        this.db = db;
        this.type = type;
        this.taxon = taxon;
        shortNames = new HashSet<String>();
        otherNames = new HashSet<String>();
    }

    public String getAccession() { return accession; }
    public int getID() { return id; }
    public String getDB() { return db; }
    public String getType() { return type; }
    public int getTaxon() { return taxon; }
    
    public Collection<String> getShortNames() { return shortNames; }
    public Collection<String> getOtherNames() { return otherNames; }
    
    public boolean containsName(String n) { 
        return accession.equals(n) || shortNames.contains(n) || otherNames.contains(n);
    }
    
    public boolean findName(String n) {
        if(accession.indexOf(n) != -1) { return true; }
        for(String sn : shortNames) { if(sn.indexOf(n) != -1) { return true; } }
        for(String on : otherNames) { if(on.indexOf(n) != -1) { return true; } }
        return false;
    }
    
    public void addShortName(String n) { shortNames.add(n); }
    public void addOtherName(String n) { otherNames.add(n); }
    
    public String getKey() { 
        return db + "/" + id + "/" + accession;
    }
    
    public int hashCode() { 
        int code = 17;
        code += id; code *= 37;
        code += accession.hashCode(); code *= 37;
        code += db.hashCode(); code *= 37;
        code += taxon; code *= 37;
        code += type.hashCode(); code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof BindInteractor)) { return false; }
        BindInteractor bi = (BindInteractor)o;
        if(id != bi.id || taxon != bi.taxon) { return false; }
        if(!accession.equals(bi.accession)) { return false; }
        if(!db.equals(bi.db)) { return false; }
        if(!type.equals(bi.type)) { return false; }
        return true;
    }
    
    public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append(db + "/" + type + "/" + id + "/" + accession);
		for(String sn : shortNames) { sb.append("/" + sn); }
		for(String on : otherNames) { sb.append("/" + on); }
		return sb.toString();
    }
}
