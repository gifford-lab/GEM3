package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;

import java.util.*;

public class RegionsWithNearbyEventPairs implements Filter<Region,Region>, Comparator<Region> {

    private Expander<Region,Region> first, second;
    private int distanceLimit;

    public RegionsWithNearbyEventPairs(Expander<Region,Region> first,
                                       Expander<Region,Region> second,
                                       int distanceLimit) {
        this.first = first;
        this.second = second;
        this.distanceLimit = distanceLimit;
    }

    public Region execute(Region region) {
        Iterator<Region> i1, i2;
        Region r1, r2;
        ArrayList<Region> al1, al2;
        Region[] a1, a2;
        al1 = new ArrayList<Region>();
        al2 = new ArrayList<Region>();
        i1 = first.execute(region);
        i2 = second.execute(region);
        while (i1.hasNext()) {
            al1.add(i1.next());
        }
        while (i2.hasNext()) {
            al2.add(i2.next());
        }
        a1 = al1.toArray(new Region[0]);
        a2 = al2.toArray(a1);
        Arrays.sort(a1,this);
        Arrays.sort(a2,this);

        int i = 0, j = 0;
        while (i < a1.length && j < a2.length ) {
            r1 = a1[i];
            r2 = a2[j];
            if (overlap(r1,r2)) {
                return region;
            }
            if (r1.getStart() < r2.getStart()) {
                i++;
            } else {
                j++;
            }
        }
        return null;
    }

    /* assume regions are on the same chromosome and genome
       and return true if they overlap */
    private boolean overlap(Region r1, Region r2) {
        return ((r1.getStart() >= (r2.getStart() - distanceLimit) && 
                 r1.getStart() <= (r2.getEnd() + distanceLimit)) ||
                (r1.getEnd() >= (r2.getStart() - distanceLimit) && 
                 r1.getEnd() <= (r2.getEnd() + distanceLimit)) ||
                (r2.getStart() >= (r1.getStart() - distanceLimit) && 
                 r2.getStart() <= (r1.getEnd() + distanceLimit)) ||
                (r2.getEnd() >= (r1.getStart() -distanceLimit) && 
                 r2.getEnd() <= (r1.getEnd()+distanceLimit)));
    }

    public int compare(Region r1, Region r2) {
        return (r1.getStart() <= r2.getStart()) ? 
            -1 : 1;
    }


}
