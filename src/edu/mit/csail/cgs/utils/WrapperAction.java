package edu.mit.csail.cgs.utils;

import java.lang.*;
import java.util.*;
import java.io.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class WrapperAction 
    extends AbstractAction { 

    private Action fBase;

    public WrapperAction(String name, Action base) { 
	super(name);
	fBase = base;
    }

    public void actionPerformed(ActionEvent e) { 
	fBase.actionPerformed(e);
    }
}
