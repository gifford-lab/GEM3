package edu.mit.csail.cgs.echo.gui;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import edu.mit.csail.cgs.echo.*;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.utils.*;

public class SelfDescribingVerbFactory implements Factory<SelfDescribingVerb> {
    
    private static final Class[] stringTypes = { String.class };

	private Class verbClass;

	public SelfDescribingVerbFactory(Class c) { 
		if(!Reverb.isSubclass(c, SelfDescribingVerb.class)) { 
			throw new IllegalArgumentException(c.getName());
		}
		verbClass = c;
	}
    
    public SelfDescribingVerb createObject(String p) { 
        
        /*
         * This is a little hackish and complicated, but the basic idea is simple.
         * If the p argument is either (null) or an empty string, then we try to 
         * invoke the nullary constructor.  If the p argument is not null, 
         * then either
         *  (a) it doesn't have a constructor which takes a string argument, and 
         *      so we just invoke the nullary constructor, or
         *  (b) it has a constructor which takes a string, and so we invoke *that* 
         *      constructor with the string p as an argument.
         *      
         * If any exception occurs, we return null.
         */
        
        try { 
            SelfDescribingVerb sdv = null;
            if(p == null || p.trim().length() == 0) { 
                sdv = (SelfDescribingVerb)verbClass.newInstance();
            } else { 
                
                Constructor cstr = null;
                try {
                    cstr = verbClass.getConstructor(stringTypes);
                } catch (SecurityException e) {
                    e.printStackTrace();
                } catch (NoSuchMethodException e) {
                    // do nothing.
                }

                if(cstr == null) { 
                    sdv = (SelfDescribingVerb)verbClass.newInstance();                    
                } else { 
                    Object[] args = new Object[1];
                    args[0] = p;
                    sdv = (SelfDescribingVerb)cstr.newInstance(args);
                }
            }
            
            return sdv;
            
        } catch(InstantiationException e) { 
            e.printStackTrace(System.err);
        } catch(IllegalAccessException e) { 
            e.printStackTrace(System.err);
        } catch (InvocationTargetException e) {
            e.printStackTrace();
        }        
        return null;
    }

	public SelfDescribingVerb createObject() {
        return createObject("");
	}
}
