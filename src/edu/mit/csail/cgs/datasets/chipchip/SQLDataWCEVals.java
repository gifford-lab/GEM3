package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;

public class SQLDataWCEVals extends SQLData {

    public SQLDataWCEVals (String exptname, String exptversion, int genomeid, Set<String> replicates) throws NotFoundException {
        super(exptname,exptversion,genomeid,replicates);
    }

    public double getMax() {
        return 65535;
    }
    public double getMax(String chrom, int start, int stop) throws NotFoundException {
        return 65535;
    }

    public double getValue(int i, int j) {
        return getWCE(i,j);
    }

}
