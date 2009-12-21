package edu.mit.csail.cgs.warpdrive;

/**
 * This class can be used to signal error cases in the Warp Drive program and 
 * as a wrapper for exceptions that may get generated within the program  
 * @author Bob
 *
 */
public class WarpDriveException extends Exception {

	/**
	   * Constructs a <code>WarpDriveException</code> with no detail message.
	   */
	  public WarpDriveException() {
	    super();
	  }

	  /**
	   * Constructs a <code>WarpDriveException</code> with specified message.
	   * @param   msg   the detail message.
	   */
	  public WarpDriveException(String msg) {
	    super(msg);
	  }

	  /**
	   * Constructs a <code>WarpDriveException</code> with an Exception Object.
	   * @param   e   the Exception Object.
	   */
	  public WarpDriveException(Exception e) {
	    super(e.getMessage());
	  }
}
