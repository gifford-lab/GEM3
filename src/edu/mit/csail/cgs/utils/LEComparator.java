package edu.mit.csail.cgs.utils;

import java.util.Comparator;

public class LEComparator implements Comparator {
    public int compare(Object o1, Object o2) {
        if (((Double)o1) <= ((Double)o2)) {
            return 1;
        } else {
            return -1;
        }
    }
    public boolean equals(Object o1) {
        return false;
    }
}

