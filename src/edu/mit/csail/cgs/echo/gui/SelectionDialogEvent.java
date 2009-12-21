/*
 * Created on Feb 21, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.gui;

import java.util.*;

public class SelectionDialogEvent extends EventObject {
    
    private SelectionDialog selDlg;

    public SelectionDialogEvent(SelectionDialog dlg) { 
        super(dlg);
        selDlg = dlg;
    }
    
    public SelectionDialog getDialog() { return selDlg; }
}
