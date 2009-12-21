package edu.mit.csail.cgs.ewok.verbs;

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseException;

/**
 *Returns the set of tiled regions in a ChipChipData object.  This is
 * similar to TiledRegionGenerator, except that by working on a
 * ChipChipData, the results may include regions tiled by multiple
 * array designs.  This is useful for determining the set of regions
 * for which results exist in a particular chipchip dataset 
*/

public class DatasetTiledRegionGenerator implements Expander<Region,Region> {

    private ChipChipData data;
    private int spacing, mincount;

    public DatasetTiledRegionGenerator (ChipChipData data, int spacing, int mincount) {
        this.data = data;
        this.spacing = spacing;
        this.mincount = mincount;
    }

    /* throws DatabaseException if the ChipChipData can't find the specified
       chromosome */
    public Iterator<Region> execute(Region r)  {
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
        } catch (NotFoundException e) {
            throw new DatabaseException(e.toString(),e);
        }


        ArrayList<Region> output = new ArrayList<Region>();
        int count = 0;
        int start = -2 * spacing;
        int last = start;
        for (int i = 0; i < data.getCount(); i++) {
            if (data.getPos(i) > last + spacing) {
                if (count >= mincount &&
                    last > 0) {
                    output.add(new Region(r.getGenome(),
                                          r.getChrom(),
                                          start,
                                          last));
                }
                count = 1;
                last = start = data.getPos(i);                
            } else {
                count++;
                last = data.getPos(i);
            }
        }
        return output.iterator();
    }
}