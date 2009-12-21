package edu.mit.csail.cgs.tools.function;

/**
 * Filters an input stream by category.  Each line should contain
 * one string.  That string is looked up in the functional annotations
 * and is present in the output iff it is assigned to one of the categories
 * specified on the command line.  If --printall is given, then
 * the program instead prints all categories retrieved for the specified gene
 *
 * java FunctionGrep --version "mus musculus" --category "Transcription Factor"
 * java FunctionGrep --version "mus musculus" --printall
 *
 * For the GO function loader (used here), the version is "genus species"
 */

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.tools.utils.Args;
import java.io.*;
import java.util.*;

public class FunctionGrep {

    public static void main(String args[]) throws Exception {
        String version = Args.parseString(args,"version",null);
        Collection<String> categories = Args.parseStrings(args,"category");
        boolean printall = Args.parseFlags(args).contains("printall");

        FunctionLoader loader = new GOFunctionLoader(GOFunctionLoader.getDefaultDBName());
        FunctionalUtils func = new FunctionalUtils(loader, version);
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String term = null;
        HashSet<String> terms = new HashSet<String>();
        while ((term = reader.readLine()) != null) {
            terms.add(term);
        }
        reader.close();
        Map<String,Set<String>> results = func.getTotalCategories(terms);
        Set<String> output = new HashSet<String>();
        if (printall) {
            for (String t : results.keySet()) {
                System.out.println(t + "\t" + results.get(t));
            }
        } else {
            for (String c : categories) {
                if (results.containsKey(c)) {
                    output.addAll(results.get(c));
                }
            }
            for (String g : output) {
                System.out.println(g);
            }

        }
        func.close();
        loader.close();
    }

}