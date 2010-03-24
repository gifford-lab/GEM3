package edu.mit.csail.cgs.tools.function;

/**
 * Filters an input stream by category.  Each line should contain
 * one string.  That string is looked up in the functional annotations
 * and is present in the output iff it is assigned to one of the categories
 * specified on the command line.  
 *
 * java FunctionGrep --version "mus musculus" --category "Transcription Factor"
 *
 * For the GO function loader (used here), the version is "genus species"
 */

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.tools.utils.Args;
import java.io.*;
import java.util.*;
import java.sql.SQLException;

public class FunctionGrep {

    public static void main(String args[]) throws Exception {
        String version = Args.parseString(args,"version",null);
        Collection<String> categories = Args.parseStrings(args,"category");
        boolean recurse = Args.parseFlags(args).contains("recurse");

        GOFunctionLoader loader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionVersion fv = loader.getVersion(version);
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String term = null;

        Set<String> acceptedNames  = new HashSet<String>();
        for (String c : categories) {
            try {
                Category category = loader.getCategory(fv,c);
                if (category != null) {
                    for (Assignment a : loader.getAssignments(category,fv)) {
                        acceptedNames.add(a.getObject());
                    }
                    if (recurse) {
                        for (Category subcat : loader.getChildCategories(category)) {
                            for (Assignment a : loader.getAssignments(subcat,fv)) {
                                acceptedNames.add(a.getObject());
                            }                        
                        }
                    }

                }
            } catch (SQLException e) {
                System.err.println("Exception trying to find " + c);
                throw e;
            }
        }
        System.err.println("Will accept " + acceptedNames);

        while ((term = reader.readLine()) != null) {
            if (acceptedNames.contains(term)) {
                System.out.println(term);
                acceptedNames.remove(term);
            }
        }
        reader.close();
        loader.close();
    }

}