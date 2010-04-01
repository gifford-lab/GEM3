package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import java.util.Map;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.io.parsing.textfiles.*;

/* fields of insert are block, col, row, name, id, type, seq */

public class AddNDFHandler extends AddGALHandler {

    public AddNDFHandler (Connection c, PreparedStatement i, Map<String,String> probeseqs) {
        super(c,i, probeseqs);
    }
    public String getIDLabel() {return "PROBE_ID";}
    public String getSequenceLabel() {return "PROBE_SEQUENCE";}
    public String getBlockLabel() {return null;}
    public String getColLabel() {return "X";}
    public String getRowLabel() {return "Y";}
    public String getNameLabel() {return null;}
    public String getControlTypeLabel() {return null;}

}
