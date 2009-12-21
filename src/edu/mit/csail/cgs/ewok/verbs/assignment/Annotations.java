package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.Collection;
import java.util.Vector;
import java.util.List;

import edu.mit.csail.cgs.ewok.verbs.Expander;

public interface Annotations<Item, Annotation> {

    /* returns the number of items that have been annotated */
	public int getNumItems();
	public Item getItem(int i);
    /* returns true iff the item has any annotations */
	public boolean isAnnotated(Item target);
    /* returns true iff the item has any annotations provided by the annotator specified by key  */
	public boolean isAnnotated(Item target, String key);
    /* returns a bit vector (as a string) showing whether the item is annotated by each
       annotator.  The order of the annotators is that specified by getAnnotationKeyOrder
    */
	public String getAnnotationBitVector(Item target);
    /* returns all annotations for this item */
	public Vector<Annotation> getAnnotations(Item target);
    /* returns all annotations for this item that were provided by the specified annotator */
	public Vector<Annotation> getAnnotations(Item target, String key);
	public List<String> getAnnotationKeyOrder();
    /* adds new items to the set of annotated items  */
	public void addItems(Collection<Item> newItems);
    /* adds a new annotator to this set of annotations.  Generates annotations for
       all previously added items
    */
	public void addAnnotations(String key, Expander<Item, ? extends Annotation> annotator);

}