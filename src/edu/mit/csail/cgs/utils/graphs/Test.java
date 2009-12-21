package edu.mit.csail.cgs.utils.graphs;

import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableModel;

public class Test extends JFrame { 

	public static void main(String[] args) {
		new Test();
	}
	
	private DirectedAlgorithms algorithms;
	private DirectedGraph graph;
	private GraphTableModel tableModel;
	
	public Test() { 
		super("Graph Tester");
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		algorithms = null;
		graph = null;
		
		tableModel = new GraphTableModel();
		JTable table = new JTable(tableModel);
		
		JPanel tablePanel = new JPanel();
		tablePanel.setLayout(new BorderLayout());
		tablePanel.add(new JScrollPane(table));
		tablePanel.setBorder(new TitledBorder("Graph"));
		
		c.add(tablePanel, BorderLayout.CENTER);
		
		setJMenuBar(createMenuBar());
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setVisible(true);
		pack();
	}
	
	private JMenuBar createMenuBar() { 
		JMenuBar bar = new JMenuBar();
		
		JMenu menu;
		JMenuItem item;
		
		bar.add(menu = new JMenu("File"));
		menu.add(new JMenuItem(getLoadGraphAction()));
		menu.add(item = new JMenuItem("Exit"));
		item.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				System.exit(0);
			}
		});
		
		bar.add(menu = new JMenu("View"));
		menu.add(item = new JMenuItem(getVisualizerAction()));
		
		bar.add(menu = new JMenu("Algorithms"));
		menu.add(new JMenuItem(getCycleCheckingAction()));
		menu.add(new JMenuItem(getTopologicalOrderingAction()));
		
		return bar;
	}
	
	public Action getLoadGraphAction() { 
		return new AbstractAction("Open...") { 
			public void actionPerformed(ActionEvent e) { 
				loadGraphFromFile();
			}
		};
	}
	
	public Action getVisualizerAction() { 
		return new AbstractAction("Visualize Graph") { 
			public void actionPerformed(ActionEvent e) { 
				GraphVisualizer viz = new GraphVisualizer(graph);
				viz.putInFrame();
			}
		};
	}
	
	public Action getCycleCheckingAction() { 
		return new AbstractAction("Check for Cycle") { 
			public void actionPerformed(ActionEvent e) { 
				if(graph != null) { 
					CycleChecker checker = new DirectedCycleChecker(graph);
					String message = checker.containsCycle() ? 
							"Graph contains cycle" : "Graph is acyclic";
					JOptionPane.showMessageDialog(Test.this.getContentPane(), 
							message, "Cycle Check", JOptionPane.INFORMATION_MESSAGE);
				}
			}
		};
	}
	
	public Action getTopologicalOrderingAction() { 
		return new AbstractAction("Topological Ordering") { 
			public void actionPerformed(ActionEvent e) { 
				if(graph != null) { 
					CycleChecker checker = new DirectedCycleChecker(graph);
					String message = null;
					String title = null;
					if(checker.containsCycle()) {
						title = "Error";
						message = "Graph contains cycle";
					} else { 
						Vector<String> ordering = algorithms.getTopologicalOrdering();
						StringBuilder sb = new StringBuilder();
						for(int i = 0; i < ordering.size(); i++) { 
							sb.append((i > 0 ? "," : "") + ordering.get(i));
						}
						title = "Topological Ordering";
						message = sb.toString();
					}
					JOptionPane.showMessageDialog(Test.this.getContentPane(), 
							message, title, JOptionPane.INFORMATION_MESSAGE);
				}
			}
		};		
	}
	
	public void loadGraphFromFile() { 
		JFileChooser chooser = new JFileChooser();
		int retval = chooser.showOpenDialog(this);
		if(retval == JFileChooser.APPROVE_OPTION) { 
			File f = chooser.getSelectedFile();
			Parser p = new Parser(f);
			try {
				DirectedGraph dg = p.parseDirectedGraph();
				setGraph(dg);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void setGraph(DirectedGraph dg) {
		graph = dg;
		algorithms = new DirectedAlgorithms(graph);
		tableModel.setGraph(dg);
	}
}

class GraphTableModel implements TableModel {
	
	private LinkedList<TableModelListener> listeners;
	private Vector<String> vertices;
	private DirectedGraph graph;
	
	public GraphTableModel() { 
		vertices = new Vector<String>();
		listeners = new LinkedList<TableModelListener>();
		graph = null;
	}
	
	public void setGraph(DirectedGraph dg) { 
		graph = dg;
		vertices = new Vector<String>(new TreeSet<String>(dg.getVertices()));
		TableModelEvent evt = new TableModelEvent(this);
		for(TableModelListener tml : listeners) { 
			tml.tableChanged(evt);
		}
	}

	public void addTableModelListener(TableModelListener tml) {
		listeners.addLast(tml);
	}

	public Class<?> getColumnClass(int c) {
		switch(c) { 
		case 0: return String.class;
		case 1: return String.class;
		default: return null;
		}
	}

	public int getColumnCount() {
		return 2;
	}

	public String getColumnName(int c) {
		switch(c) { 
		case 0: return "Vertex";
		case 1: return "Parents";
		default: return null;
		}
	}

	public int getRowCount() {
		return vertices.size();
	}

	public Object getValueAt(int r, int c) {
		switch(c) { 
		case 0: return vertices.get(r);
		case 1: 
			StringBuilder sb = new StringBuilder();
			if(graph != null) { 
				for(String parent : graph.getParents(vertices.get(r))) { 
					sb.append(parent + ",");
				}
			}
			return sb.toString();
		default:
			return null;
		}
	}

	public boolean isCellEditable(int arg0, int arg1) {
		return false;
	}

	public void removeTableModelListener(TableModelListener tml) {
		listeners.remove(tml);
	}

	public void setValueAt(Object arg0, int arg1, int arg2) {
	} 
	
}