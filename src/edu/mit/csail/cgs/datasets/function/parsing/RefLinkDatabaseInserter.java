/*
 * Created on Feb 13, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.io.*;
import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class RefLinkDatabaseInserter {
    
    private Genome genome;

    public RefLinkDatabaseInserter(Genome g) { 
        genome = g;
    }
    
    public void insertIntoDB(String versionName) throws SQLException, UnknownRoleException {
        DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
        FunctionVersion version = loader.getVersion(versionName);
        
        java.sql.Connection cxn = genome.getUcscConnection();
        Statement s = cxn.createStatement();
        
        Map<String,Set<String>> assigns = new HashMap<String,Set<String>>();
        
        System.out.println("Building Assignment Map...");
        ResultSet rs = s.executeQuery("select r.protAcc, r.mrnaAcc, r.locusLinkId, k.spID, k.kgID " +
                "from refLink r, kgXref k where r.protAcc=k.protAcc");
        int i = 0;
        while(rs.next()) { 
            String prot = rs.getString(1);
            String mrna = rs.getString(2);
            String entrez = rs.getString(3);
            String spID = rs.getString(4);
            String kgID = rs.getString(5);
            
            Collection<Assignment> protAssigns = loader.getAssignments(spID, version);
            for(Assignment a : protAssigns) { 
                Category cat = a.getCategory();
                String goID = cat.getName();
                if(!assigns.containsKey(goID)) { assigns.put(goID, new HashSet<String>()); }
                
                assigns.get(goID).add(prot);
                assigns.get(goID).add(mrna);
                assigns.get(goID).add(entrez);
                assigns.get(goID).add(kgID);
            }
            
            i++;
            if(i % 1000 == 0) { System.out.print("."); System.out.flush(); }
            if(i % 10000 == 0) { System.out.print("[" + i + "]"); System.out.flush(); }
        }
        System.out.println();
        
        rs.close();
        s.close();
        
        System.out.println("Inserting Assignments...");
        //loader.insertAssignments(versionName, assigns);
        
        DatabaseFactory.freeConnection(cxn);
        loader.close();
    }
}
