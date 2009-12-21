/**
 * 
 */
package edu.mit.csail.cgs.utils.graphs.ui;

import java.util.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import edu.mit.csail.cgs.utils.graphs.*;

/**
 * @author Timothy Danford
 *
 */
public class GraphEditor extends JPanel {
	
	public static void main(String[] args) { 
		DirectedGraph dg = new DirectedGraph();
		dg.addVertex("node1");
		dg.addVertex("node2");
		dg.addVertex("node3");
		dg.addEdge("node1", "node2");
		dg.addEdge("node1", "node3");
		dg.addEdge("node2", "node3");
		
		Frame f = new Frame(dg);
	}
	
	public static class Frame extends JFrame {
		private Graph graph;
		
		public Frame(Graph g) { 
			super("Graph Editor");
			graph = g;
			
			GraphEditor ge = new GraphEditor(graph);
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(ge, BorderLayout.CENTER);
			
			setVisible(true);
			pack();
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
		
	}

	private JList nodeList, edgeList;
	private DefaultListModel nodeModel, edgeModel;
	private Graph graph;
	private boolean directed;
	private JTextField entryField;
	private JButton addVertex, addEdge, deleteEdge, deleteVertex, viewGraph;
	
	public GraphEditor(Graph g) { 
		graph = g;
		directed = graph instanceof DirectedGraph;
		nodeModel = new DefaultListModel();
		edgeModel = new DefaultListModel();
		nodeList = new JList(nodeModel);
		edgeList = new JList(edgeModel);
		
		JPanel entryPanel = new JPanel(); entryPanel.setLayout(new GridLayout(1, 6));
		entryPanel.add(entryField = new JTextField());
		entryPanel.add(addVertex = new JButton("+V"));
		entryPanel.add(deleteVertex = new JButton("-V"));
		entryPanel.add(addEdge = new JButton("+E"));
		entryPanel.add(deleteEdge = new JButton("-E"));
		entryPanel.add(viewGraph = new JButton("View"));
		
		setLayout(new BorderLayout());
		JPanel viewPanel = new JPanel(); viewPanel.setLayout(new GridLayout(1, 2));
		viewPanel.add(new JScrollPane(nodeList));
		viewPanel.add(new JScrollPane(edgeList));
		
		add(viewPanel, BorderLayout.CENTER);
		add(entryPanel, BorderLayout.SOUTH);
		
		viewGraph.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				JFrame f = new Visualizer.Frame(graph);
				f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			}
		});
		
		addVertex.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) {
				String name = entryField.getText();
				if(!graph.getVertices().contains(name)) { 
					if(directed) { 
						DirectedGraph dg = (DirectedGraph)graph;
						dg.addVertex(name);
					} else { 
						UndirectedGraph ug = (UndirectedGraph)graph;
						ug.addVertex(name);
					}
					nodeModel.addElement(name);
				}
			}
		});
		deleteVertex.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				int[] selected = nodeList.getSelectedIndices();
				if(selected.length == 1) {
					String node = (String)nodeModel.get(selected[0]);
					if(graph.getVertices().contains(node)) { 
						if(directed) { 
							DirectedGraph dg = (DirectedGraph)graph;
							dg.removeVertex(node);
						} else { 
							UndirectedGraph ug = (UndirectedGraph)graph;
							ug.removeVertex(node);
						}
						nodeModel.remove(selected[0]);
						edgeModel.clear();
					}
				}
			}
		});
		addEdge.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				int[] selected = nodeList.getSelectedIndices();
				if(selected.length == 1) { 
					String start = (String)nodeModel.get(selected[0]);
					String target = entryField.getText();
					if(graph.getVertices().contains(target) && 
							graph.getVertices().contains(start) && 
							!graph.getNeighbors(start).contains(target)) { 
						if(directed) { 
							DirectedGraph dg = (DirectedGraph)graph;
							dg.addEdge(start, target);
						} else { 
							UndirectedGraph ug = (UndirectedGraph)graph;
							ug.addEdge(start, target);
						}
						edgeModel.addElement(target);
					}
				}
			}
		});
		deleteEdge.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) {
				int[] selNodes = nodeList.getSelectedIndices();
				if(selNodes.length == 1) { 
					String node = (String)nodeModel.get(selNodes[0]);
					int[] selEdges = edgeList.getSelectedIndices();
					for(int i = selEdges.length-1; i >= 0; i--) { 
						String target = (String)edgeModel.get(selEdges[i]);
						if(directed) { 
							DirectedGraph dg = (DirectedGraph)graph;
							dg.removeEdge(node, target);
						} else { 
							UndirectedGraph ug = (UndirectedGraph)graph;
							ug.removeEdge(node, target);
						}
						edgeModel.remove(selEdges[i]);
					}
				}
			}
		});
		
		for(String v : new TreeSet<String>(graph.getVertices())) {  
			nodeModel.addElement(v);
		}

		nodeList.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				edgeModel.clear();
				int[] selected = nodeList.getSelectedIndices();
				if(selected.length == 1) { 
					String value = (String)nodeModel.get(selected[0]);
					TreeSet<String> neighbors = new TreeSet<String>(graph.getNeighbors(value));
					for(String n : neighbors) { 
						edgeModel.addElement(n);
					}
				}
			} 
		});
	}
}
