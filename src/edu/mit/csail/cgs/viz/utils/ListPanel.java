/*
 * Created on Aug 20, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.Collection;

/**
 * @author tdanford
 */
public interface ListPanel<X> extends ListPanelEventSource {
    public int getNumValues();
    public X getValue(int index);
    public void addValue(X v);
    public void clear();
    public X getFirstSelectedValue();
    public Collection<X> getSelectedValues();
    public void removeValue(int index);
    public int[] getSelectedIndices();
    public void removeSelectedValues();
    public Collection<X> getAllValues();
}