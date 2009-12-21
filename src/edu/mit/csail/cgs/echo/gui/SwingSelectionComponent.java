package edu.mit.csail.cgs.echo.gui;

import javax.swing.JComponent;

public interface SwingSelectionComponent<X> {
	public X getSelectedValue();
	public JComponent asJComponent();
}
