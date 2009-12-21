package edu.mit.csail.cgs.tools.core;

import java.sql.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.database.Sequence;

/**
 * adds the listed factors to the DB:
 * AddFactors H3K4me3 'polyA-RNA'
 */

public class AddFactors {

    public static void main(String args[]) throws Exception {
        java.sql.Connection cxn = DatabaseFactory.getConnection("core");        
        PreparedStatement insert = cxn.prepareStatement("insert into factors(id,name) values(" +
                                                        Sequence.getInsertSQL(cxn,"factors_id") + ",?)");
        for (int i = 0; i < args.length; i++) {
            insert.setString(1,args[i]);
            try {
                insert.execute();
            } catch (SQLException e) {
                System.err.println(e.toString());
            }

        }

    }

}