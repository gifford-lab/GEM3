/*
 * Created on May 24, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;
import edu.mit.csail.cgs.utils.Pair;

/**
 * @author tdanford
 */
public interface FunctionalCategories {
    public Collection<String> getCategories();
    public String getDescription(String name);
    public Collection<String> getParents(String name);
    public Collection<String> classify(String object);
	public Collection<String> getCategoryObjects(String catName);
	public boolean hasObject(String obj);
	public boolean hasCategory(String cat);
}

