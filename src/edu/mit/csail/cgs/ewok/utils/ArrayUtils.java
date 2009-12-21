package edu.mit.csail.cgs.ewok.utils;

import java.util.*;

public class ArrayUtils {

    public static <T> ArrayList<T> iteratorToList(Iterator<T> i) {
        ArrayList<T> list = new ArrayList<T>();
        while (i.hasNext()) {
            list.add(i.next());
        }
        return list;
    }

}
