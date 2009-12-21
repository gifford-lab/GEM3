package edu.mit.csail.cgs.tools.function;

/**
 * Returns a list of all categories in a set of annotations that
 * match a regular expression.
 *
 * java CategoryGrep --version "mus musculus" --regex "transcription"
 */

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.tools.utils.Args;
import java.util.*;
import java.util.regex.*;


public class CategoryGrep {

    public static void main(String args[]) throws Exception {
        String version = Args.parseString(args,"version",null);
        String regex = Args.parseString(args,"regex","");
        FunctionLoader loader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionVersion fv = loader.getVersion(version);
        Collection<Category> cats = loader.getCategories(fv);
        Pattern patt = Pattern.compile(regex);
        for (Category c : cats) {
            String n = c.getName();
            Matcher m = patt.matcher(n);
            if (m.find()) {
                System.out.println(n);
            }
        }
        loader.close();
    }

        


}
