/*
 * Created on Apr 1, 2006
 */
package edu.mit.csail.cgs.ewok.utils;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.orthology.*;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class EwokOrthology extends EwokBase {
    
    private OrthologyLoader loader;
    private OrthologyMapping mapping;

    public EwokOrthology(OrthologyLoader loader, String mapName, String mapVersion) throws NotFoundException {
        super();
        this.loader = loader;
        mapping = loader.loadMapping(mapName, mapVersion);
    }

    public EwokOrthology(String sp, String gn, OrthologyLoader loader, String mapName, String mapVersion) throws NotFoundException {
        super(sp, gn);
        this.loader = loader;
        mapping = loader.loadMapping(mapName, mapVersion);
    }
    
    public Iterator<String> translateNames(Genome otherGenome, Iterator<String> baseNames) { 
        OrthologyExpander exp = new OrthologyExpander(loader, mapping, getGenome(), otherGenome);
        Iterator<String> secondNames = new ExpanderIterator<String,String>(exp, baseNames);
        return secondNames;
    }
    
    public static class StringsFromInputStream implements Iterator<String> { 
        
        private BufferedReader reader;
        private String nextLine; 
        
        public StringsFromInputStream(InputStream is) throws IOException { 
            reader = new BufferedReader(new InputStreamReader(is));
            nextLine = reader.readLine();
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            return nextLine != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public String next() {
            String ret = nextLine;
            nextLine = null;
            try {
                nextLine = reader.readLine();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return ret;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public static class GeneNameExpander implements Expander<Gene,String> {
        
        public GeneNameExpander() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(java.lang.Object)
         */
        public Iterator<String> execute(Gene a) {
            LinkedList<String> names = new LinkedList<String>();
            names.addLast(a.getName());
            names.addLast(a.getID());
            return names.iterator();
        } 
        
    }
}
