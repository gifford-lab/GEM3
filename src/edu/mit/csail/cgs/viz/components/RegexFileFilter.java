package edu.mit.csail.cgs.viz.components;

import java.io.File;
import javax.swing.filechooser.FileFilter;

public class RegexFileFilter extends FileFilter {

    private String regex, desc;
    private boolean dir;
    /**
     * Creates a RegexFileFilter that accepts 
     * all files whose name (but not path) matches
     * the regularExpression.
     * if acceptDirectories is true, then all directories
     * are accepted (necessary if used with JFileChooser
     * and you want to be able to navigate into directories)
     */
    public RegexFileFilter(String regularExpression,
                           String description,
                           boolean acceptDirectories) {
        desc = description;
        regex = regularExpression;
        dir = acceptDirectories;
    }
    public String getDescription() {return desc;}
    public boolean accept(File f) {
        return (dir && f.isDirectory()) || f.getName().matches(regex);
    }

}