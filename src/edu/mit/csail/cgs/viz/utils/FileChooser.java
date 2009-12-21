/*
 * Created on Aug 21, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.io.File;
import edu.mit.csail.cgs.utils.ObjectChooser;
import javax.swing.*;

/**
 * @author tdanford
 */
public class FileChooser implements ObjectChooser<File> {

    private JFrame parent;
    private String name;
    private JFileChooser chooser;

    public FileChooser(JFrame p) {
        parent = p;
        name = null;
        chooser = new JFileChooser();
    }
    
    public FileChooser(JFrame p, String n) {
        parent = p;
        name = n;
        chooser = new JFileChooser();
        chooser.setApproveButtonText(n);
    }
    
    public File choose() { 
    	return chooseOpen();
    }

    public File chooseSave() {
        chooser.setMultiSelectionEnabled(false);
        int option = chooser.showSaveDialog(parent);
        if(option == JFileChooser.APPROVE_OPTION && chooser.getSelectedFile() != null) {
            return chooser.getSelectedFile();
        } else { 
            return null;
        }
    }

    public File chooseOpen() {
        chooser.setMultiSelectionEnabled(false);
        int option = chooser.showOpenDialog(parent);
        if(option == JFileChooser.APPROVE_OPTION && chooser.getSelectedFile() != null) {
            return chooser.getSelectedFile();
        } else { 
            return null;
        }
    }

    public File[] chooseAll() {
        chooser.setMultiSelectionEnabled(true);
        int option = chooser.showOpenDialog(parent);
        if(option == JFileChooser.APPROVE_OPTION && chooser.getSelectedFile() != null) {
            return chooser.getSelectedFiles();
        } else { 
            return null;
        }
    }

}
