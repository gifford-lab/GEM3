package edu.mit.csail.cgs.tools.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;

/**
 * Dumps data from the alignmentparameters table for a single alignment 
 * 
 * GetAlignmentParameters --id ID
 */

public class GetAlignmentParameters {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        ChipSeqLoader loader = new ChipSeqLoader();
        
        Integer id = Args.parseInteger(args,"id", -1);
        if(id != -1){
        	ChipSeqAlignment align = loader.loadAlignment(id);
        	Map<String,String> params = loader.getAlignmentParameters(align);
        	
        	for(String s : params.keySet()){
        		System.out.println(s+"="+params.get(s));
        	}
        }
    }
}