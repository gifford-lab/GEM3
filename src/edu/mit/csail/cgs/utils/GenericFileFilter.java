package edu.mit.csail.cgs.utils;

import java.io.File;
import javax.swing.filechooser.FileFilter;

/**
 * This is a file filter for use with a JFileChooser. It limits the JFileChooser
 * to only show/allow files with a specified set of extensions
 * 
 * @author Bob
 *
 */
public class GenericFileFilter extends FileFilter {

	String[] extensions;
	String description;


	/**
	 *
	 * @param extensions String[] extensions to allow
	 * @param description String a description of the types of files with these extensions
	 */
	public GenericFileFilter(String[] extensions, String description) {
		this.extensions = extensions;
		this.description = description;
	}


	/**
	 *
	 * @param f File a file to test for acceptance by this filter
	 * @return boolean whether this filter accepts the specified file
	 */
	public boolean accept(File f) {
		if (f.isDirectory()) {
			return true;
		}

		String extension = getFileExtension(f.getName());
		if (extension != null) {
			boolean valid = false;
			for (int i = 0; i < extensions.length; i++) {
				if (extensions[i].equals(extension)) {
					valid = true;
					break;
				}
			}
			return valid;
		}
		else {
			return false;
		}
	}


	/**
	 *
	 * @return String
	 */
	public String getDescription() {
		return description;
	}


	/**
	 *
	 * @param filename String a name of a file
	 * @return String the extension of that file
	 */
	private String getFileExtension(String filename) {
		String extension = null;
		int i = filename.lastIndexOf('.');

		if (i > 0 && i < filename.length() - 1) {
			extension = filename.substring(i + 1).toLowerCase();
		}
		return extension;
	}
}
