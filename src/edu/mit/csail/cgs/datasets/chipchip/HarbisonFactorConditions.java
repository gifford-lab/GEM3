package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class HarbisonFactorConditions {

    public static void main(String[] args) { 
        try {
            HarbisonFactorConditions hfc = new HarbisonFactorConditions();

            for(String f : new TreeSet<String>(hfc.getFactors())) {
                System.out.println(f);
                for(String c : new TreeSet<String>(hfc.getConditions(f))) { 
                    System.out.println("\t" + c);
                }
            }

        } catch (SQLException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    private Genome yeastGenome;
    private Map<String,Set<String>> factor2Conds, cond2Factors;
    private int size;

    public HarbisonFactorConditions() throws SQLException, NotFoundException {
        yeastGenome = Organism.findGenome("sacCer1");
        factor2Conds = new HashMap<String,Set<String>>();
        cond2Factors = new HashMap<String,Set<String>>();
        size = 0;
        loadFromDB();
    }

    public Set<String> getFactors() { 
        return factor2Conds.keySet();
    }

    public Set<String> getConditions() {
        return cond2Factors.keySet();
    }

    public Set<String> getFactors(String c) { return cond2Factors.get(c); }
    public Set<String> getConditions(String f) { return factor2Conds.get(f); }

    public void loadFromDB() throws SQLException {
        size = 0;
        cond2Factors.clear(); factor2Conds.clear();

        java.sql.Connection cxn = 
            yeastGenome.getUcscConnection();
        Statement s = cxn.createStatement();

        ResultSet rs = 
            s.executeQuery("select name, growthCondition from transRegCodeCondition");

        while(rs.next()) { 
            String factor = rs.getString(1);
            String condition = rs.getString(2);
            addFactorCondition(factor,condition);
        }

        rs.close();

        s.close();
        DatabaseFactory.freeConnection(cxn);
    }

    private void addFactorCondition(String f, String c) { 
        if(!factor2Conds.containsKey(f)) { factor2Conds.put(f, new HashSet<String>()); }
        if(!cond2Factors.containsKey(c)) { cond2Factors.put(c, new HashSet<String>()); }
        if(!factor2Conds.get(f).contains(c)) { size += 1; }
        factor2Conds.get(f).add(c);
        cond2Factors.get(c).add(f);
    }
}
