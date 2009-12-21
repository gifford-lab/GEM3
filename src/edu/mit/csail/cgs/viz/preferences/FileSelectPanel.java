package edu.mit.csail.cgs.viz.preferences;

import java.io.*;
import java.awt.*;
import javax.swing.*;

import edu.mit.csail.cgs.viz.utils.FileChooser;

import java.awt.event.*;

public class FileSelectPanel extends JPanel {
	
	public static void main(String[] args) { 
		TestFrame tf = new TestFrame(null);
	}

	private File file;
	private JButton chooseButton;
	private JLabel label;
	
	public FileSelectPanel(File f) { 
		super();
		file = f;
		init();
	}
		
		
	public void init() {
		Runnable r = new Runnable() {
			public void run() { 
				label = new JLabel(file == null ? "None" : file.getName());
				chooseButton = new JButton("Choose File");
				setLayout(new GridLayout(1, 2));
				add(label);
				add(chooseButton);
				
				chooseButton.addActionListener(new ActionListener() { 
					public void actionPerformed(ActionEvent e) { 
						chooseNewFile();
					}
				});
			}
		};
		EventQueue.invokeLater(r);
	}
	
	public void chooseNewFile() { 
		FileChooser chooser = new FileChooser(null);
		file = chooser.choose();
		updateLabel();
	}
	
	private void updateLabel() { 
		Runnable r = new Runnable() { 
			public void run() { 
				label.setText(file == null ? "None" : file.getName());
				repaint();
			}
		};
		EventQueue.invokeLater(r);
	}
	
	public File getFile() { return file; }
	
}

class TestFrame extends JFrame {
	private FileSelectPanel fsp; 
	public TestFrame(File f) { 
		super("File Test");
		fsp = new FileSelectPanel(f);

		Runnable r = new Runnable() {
			public void run() { 
				Container c = (Container)getContentPane();
				c.setLayout(new BorderLayout());
				c.add(fsp, BorderLayout.NORTH);
				
				addWindowListener(new WindowAdapter() { 
					public void windowClosing(WindowEvent e) { 
						System.out.println("File: " + fsp.getFile().getName());
					}
				});
				
				setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
				setVisible(true);
				pack();
			}
		};
		
		EventQueue.invokeLater(r);
	}
}