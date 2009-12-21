/*
 * Created on Nov 9, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;

/**
 * @author tdanford
 */
public class GOALine {

    public String[] array;
    private String db, dbObjectID, dbObjectSymbol;
    private String qualifier;
    private String goID;
    private String reference, evidence, with;
    private String aspect;
    private String dbObjectName;
    private String synonym;
    private String dbObjectType;
    private String taxonID;
    private String date, assignedBy;
    
    public GOALine(String line) { 
        array = line.split("\t");
        db = array[0]; 
        dbObjectID = array[1];
        dbObjectSymbol = array[2];
        qualifier = array[3];
        goID = array[4];
        reference = array[5];
        evidence = array[6];
        with = array[7];
        aspect = array[8];
        dbObjectName = array[9];
        synonym = array[10];
        dbObjectType = array[11];
        taxonID = array[12];
        date = array[13];
        assignedBy = array[14];
    }
    
    public String getDB() { return db; }
    public String getObjectID() { return dbObjectID; }
    public String getObjectSymbol() { return dbObjectSymbol; }
    public String getObjectName() { return dbObjectName; }
    public String getObjectType() { return dbObjectType; }
    public String getTaxonID() { return taxonID; }
    public String getAspect() { return aspect; }
    public String getGOID() { return goID; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof GOALine)) { return false; }
        GOALine l = (GOALine)o;
        if(!db.equals(l.db)) { return false; }
        if(!dbObjectID.equals(l.dbObjectID)) { return false; }
        if(!dbObjectSymbol.equals(l.dbObjectSymbol)) { return false; }
        if(!goID.equals(l.goID)) { return false; }
        if(!dbObjectName.equals(l.dbObjectName)) { return false; }
        if(!taxonID.equals(l.taxonID)) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += db.hashCode(); code *= 37;
        code += dbObjectID.hashCode(); code *= 37;
        code += goID.hashCode(); code *= 37;
        code += taxonID.hashCode(); code *= 37;
        return code;
    }
    
    public String toString() { return dbObjectSymbol + " (" + dbObjectID + ") --> " + goID; }
    
    public void addToMap(Map<String,Set<String>> map) { 
        if(!map.containsKey(goID)) { map.put(goID, new HashSet<String>()); }
        map.get(goID).add(dbObjectID);
        map.get(goID).add(dbObjectSymbol);
    }
}
