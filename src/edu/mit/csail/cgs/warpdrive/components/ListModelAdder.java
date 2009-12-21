package edu.mit.csail.cgs.warpdrive.components;

import javax.swing.DefaultListModel;
import java.util.Collection;

class ListModelAdder implements Runnable {
    private Collection<Object> collection;
    private Object[] array;
    private DefaultListModel listmodel;
    public ListModelAdder(Collection o,
                          DefaultListModel l) {
        collection = o;
        listmodel = l;
        array = null;
    }
    public ListModelAdder(Object[] o,DefaultListModel l){
        collection = null;
        listmodel  = l;
        array = o;
    }
        
    public void run() {
        if (collection != null) {
            synchronized (collection) {
                for (Object o : collection) {
                    listmodel.addElement(o);
                }
                collection.clear();
            };
        }
        if (array != null) {
            for (int i = 0; i < array.length; i++) {
                listmodel.addElement(array[i]);
            }
        }
    }

}
