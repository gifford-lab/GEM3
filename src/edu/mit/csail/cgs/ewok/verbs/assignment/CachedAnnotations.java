/*
 * Created on Dec 1, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

public class CachedAnnotations<Item,Annotation> implements Annotations<Item, Annotation> {

    private Vector<Item> items;
    private Map<String,Expander<Item,? extends Annotation>> annotators;
    private Vector<String> annotatorKeyOrder;    
    private Map<Item,Map<String,Vector<Annotation>>> annotations;
    
    public CachedAnnotations() {
        items = new Vector<Item>();
        annotators = new HashMap<String,Expander<Item,? extends Annotation>>();
        annotatorKeyOrder = new Vector<String>();        
        annotations = new HashMap<Item,Map<String,Vector<Annotation>>>();
    }
    
    public CachedAnnotations(Collection<Item> tgs) {
        items = new Vector<Item>(tgs);
        annotators = new HashMap<String,Expander<Item,? extends Annotation>>();
        annotatorKeyOrder = new Vector<String>();
        
        annotations = new HashMap<Item,Map<String,Vector<Annotation>>>();
        for(Item target : items) { 
            annotations.put(target, new HashMap<String,Vector<Annotation>>());
        }
    }
    
    public int getNumItems() { return items.size(); }
    public Item getItem(int i) { return items.get(i); }
	public boolean hasItem(Item itm) { return items.contains(itm); }
    
    public String getAnnotationBitVector(Item target) { 
        StringBuilder sb = new StringBuilder();
        for(String k : annotatorKeyOrder) { 
            sb.append((isAnnotated(target, k) ? "1" : "0"));
        }
        return sb.toString();
    }
    
    public boolean isAnnotated(Item target, String key) {  
        return !annotations.get(target).get(key).isEmpty();
    }
    
    public boolean isAnnotated(Item target) { 
        for(String k : annotatorKeyOrder) { 
            if(isAnnotated(target, k)) { return true; }
        }
        return false;
    }
    
    public Vector<Annotation> getAnnotations(Item target) { 
        Vector<Annotation> evts = new Vector<Annotation>();
        for(String key : annotatorKeyOrder) { 
            evts.addAll(annotations.get(target).get(key));
        }
        return evts;
    }
    
    public Vector<Annotation> getAnnotations(Item target, String key) { 
        return annotations.get(target).get(key);
    }
    
    public List<String> getAnnotationKeyOrder() { 
        return (List<String>)annotatorKeyOrder.clone();
    }
    
    public void addItems(Collection<Item> newItems) { 
        for(Item target : newItems) {
            if(items.contains(target)) { throw new IllegalArgumentException(target.toString()); }
            items.add(target);
            annotations.put(target, new HashMap<String,Vector<Annotation>>());
            
            for(String k : annotatorKeyOrder) { 
                Expander<Item,? extends Annotation> annotator = annotators.get(k);
                Vector<Annotation> annots = new Vector<Annotation>();
                annotations.get(target).put(k, annots);
                
                Iterator<? extends Annotation> itemAnnotations = annotator.execute(target);
                while(itemAnnotations.hasNext()) { 
                    Annotation a = itemAnnotations.next();
                    annots.add(a);
                }
            }
        }
    }

    public void addAnnotations(String key, Expander<Item,? extends Annotation> annotator) { 
        if(annotators.containsKey(key)) { throw new IllegalArgumentException(key); }
        annotatorKeyOrder.add(key);
        annotators.put(key, annotator);
        int nonZeroCount = 0;
        
        for(Item target : items) { 
            Map<String,Vector<Annotation>> targetAnnotations = annotations.get(target);
            
            Iterator<? extends Annotation> evts = annotator.execute(target);
            Vector<Annotation> eventVector = new Vector<Annotation>();
            while(evts.hasNext()) { 
                Annotation evt = evts.next();
                eventVector.add(evt);
            }
            
            targetAnnotations.put(key, eventVector);
            nonZeroCount += (eventVector.isEmpty() ? 0 : 1);
            markProgress();
        }
        
        System.out.println("Annotated: \"" + key + "\" with " + nonZeroCount + " items.");
    }
    
    /* called every time an item is annotated.  Subclasses can use this
       to display a measure of progress of the annotation */
    public void markProgress() {}
}
