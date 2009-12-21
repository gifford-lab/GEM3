package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.io.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.database.*;

public class HarbisonRegCodeProbeGenerator<X extends Region> implements Expander<X,HarbisonRegCodeProbe> {
    
    public static void main(String[] args) {
        try {
            Genome g = Organism.findGenome("sacCer1");
            HarbisonRegCodeProbeGenerator<NamedRegion> gen = 
                new HarbisonRegCodeProbeGenerator<NamedRegion>(g);
            Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
            Iterator<HarbisonRegCodeProbe> probes = 
                new ExpanderIterator<NamedRegion,HarbisonRegCodeProbe>(gen, chroms);
            while(probes.hasNext()) { 
                HarbisonRegCodeProbe probe = probes.next();
                System.out.println(probe.getName());
                for(String f : probe.getFactors()) { 
                    System.out.println("\t" + f);
                    for(String c : probe.getConditions(f)) {
                        int str = probe.getBindingStrength(f, c);
                        String strString = HarbisonRegCodeProbe.strengthStrings[str];
                        System.out.println("\t\t" + c + " : " + strString);
                    }
                }
                System.out.println();
            }
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
        
    }
    
    private Genome genome;
    
    public HarbisonRegCodeProbeGenerator(Genome g) {
        genome = g;
    }

    public Iterator<HarbisonRegCodeProbe> execute(X region) {
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = 
                cxn.prepareStatement("select name, chrom, chromStart, chromEnd, " +
                        "tfCount, tfList, bindVals from transRegCodeProbe where chrom = ? and " +
                "((chromStart <= ? and chromEnd >= ?) or " +
                "(chromStart >= ? and chromStart <= ?)) order by chromStart");
            
            String chr = region.getChrom();
            if (!chr.matches("^(chr|scaffold).*")) {
                chr = "chr" + chr;
            }
            ps.setString(1,chr);
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getStart());
            ps.setInt(5,region.getEnd());
            ArrayList<HarbisonRegCodeProbe> results = new ArrayList<HarbisonRegCodeProbe>();
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                String name = rs.getString(1), chrom = rs.getString(2);
                int start = rs.getInt(3), end = rs.getInt(4), tfCount = rs.getInt(5);
                
                String tfList = "", bindVals = "";
                
                if(tfCount > 0) { 
                    Blob tfListBlob = rs.getBlob(6);
                    Blob bindValsBlob = rs.getBlob(7);

                    BufferedReader br = 
                        new BufferedReader(new InputStreamReader(tfListBlob.getBinaryStream()));
                    tfList = br.readLine();
                    //System.out.println(tfList);
                    br.close();

                    br = new BufferedReader(new InputStreamReader(bindValsBlob.getBinaryStream()));
                    bindVals = br.readLine();
                    //System.out.println(bindVals);
                    br.close();
                }
                
                if(chrom.matches("^chr.*")) { 
                	chrom = chrom.substring(3, chrom.length());
                }
                results.add(new HarbisonRegCodeProbe(region.getGenome(), 
                		chrom, start, end, name, tfCount, tfList, bindVals));
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();

        } catch (SQLException ex) {
            throw new edu.mit.csail.cgs.utils.database.DatabaseException("Couldn't get UCSC RefGenes",ex);
        } catch (IOException e) {
            throw new edu.mit.csail.cgs.utils.database.DatabaseException("Error Reading the tfList/bindVals blobs", e);
        }
    }
}    
