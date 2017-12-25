package edu.mit.csail.cgs.viz.components;

import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.table.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqExpt;
import edu.mit.csail.cgs.datasets.species.Genome;

public abstract class GenericSelectPanel<X> extends JSplitPane implements Closeable, Runnable {

	private JButton addButton, removeButton, filterButton;

	protected JTable filteredList, selectedList;
	protected ObjectTableModel<X> filteredModel, selectedModel;
	protected JPanel buttonPanel;
	private Genome genome;
	protected HashMap<String, ChipSeqExpt> readdb;				// the readdb meta file matching the genome
	private boolean handlingNewGenome, dataReady;
	private Thread thread;
	/*
	 * asyncSelected is here so you can call addToSelected() while
	 * updateComponents() is running. Since updateComponents probably starts by
	 * clearing them, your changes will be lost. By sticking things in asyncSelected
	 * and then adding them at the end, your additions will be kept
	 */
	private ArrayList<X> asyncSelected;

	public GenericSelectPanel() {
		super(JSplitPane.VERTICAL_SPLIT);
		setDividerLocation(200);
		handlingNewGenome = false;
		dataReady = false;
		asyncSelected = new ArrayList<X>();
	}

	public GenericSelectPanel(Genome g) {
		super(JSplitPane.VERTICAL_SPLIT);
		setDividerLocation(200);
		handlingNewGenome = false;
		dataReady = false;
		genome = g;
		asyncSelected = new ArrayList<X>();
	}

	public Genome getGenome() {
		return genome;
	}

	/*
	 * sets the genome and starts a new thread to process that change. Database work
	 * should happen in that thread and the results should be cached. Swing
	 * components should be updated with a call to SwingUtilities.invokeLater().
	 * Since this also invokes run(), the run() method must determine the current
	 * operation by querying: handlingNewGenome() and dataReady()
	 * 
	 */
	public void setGenome(Genome g,  HashMap<String, ChipSeqExpt> readdb) {
		if (g == null) {
			return;
		}
		genome = g;
		this.readdb = readdb;
		if (thread != null && thread.isAlive()) {
			thread.interrupt();
			thread = null;
		}

		thread = new Thread(this);
		handlingNewGenome = true;
		dataReady = false;
		thread.start();
	}
	public void setGenome(Genome g) {	// just a dummy
	}

	public boolean handlingNewGenome() {
		return handlingNewGenome;
	}

	public boolean dataReady() {
		return dataReady;
	}

	public void run() {
		if (handlingNewGenome()) {
			retrieveData();
			handlingNewGenome = false;
			dataReady = true;
			SwingUtilities.invokeLater(this);
		} else if (dataReady()) {
			updateComponents();
			synchronized (asyncSelected) {
				for (X x : asyncSelected) {
					selectedModel.addObject(x);
				}
				asyncSelected.clear();
				dataReady = false;
			}
		}
	}

	/*
	 * retrieve data from the database; must cache somewhere
	 */
	public abstract void retrieveData();

	/*
	 * called from swing thread to update the swing components for the cached data
	 */
	public abstract void updateComponents();

	/* sub-class constructors must call one of the two init() methods */
	public void init(JTable filteredList, JTable selectedList, ObjectTableModel<X> filteredModel,
			ObjectTableModel<X> selectedModel) {
		final JTableHeader header = filteredList.getTableHeader();
		final ObjectTableModel filteredmodel = filteredModel;
		header.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				filteredmodel.sortByColumn(header.columnAtPoint(e.getPoint()));
			}
		});

		this.filteredList = filteredList;
		this.selectedList = selectedList;
		this.filteredModel = filteredModel;
		this.selectedModel = selectedModel;
		addButton = new JButton("Add");
		removeButton = new JButton("Remove");
		filterButton = new JButton("Filter");
		final GenericSelectPanel panel = this;
		addButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				panel.add();
			}
		});
		removeButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				panel.remove();
			}
		});
		filterButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				panel.filter();
			}
		});
		filteredList.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				panel.filterToSel(e);
			}
		});
		selectedList.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				panel.delSel(e);
			}
		});
		JPanel selPanel = new JPanel();
		selPanel.setLayout(new BorderLayout());
		JPanel actPanel = new JPanel();
		actPanel.setLayout(new BorderLayout());
		selPanel.add(new JScrollPane(selectedList), BorderLayout.CENTER);
		actPanel.add(new JScrollPane(filteredList), BorderLayout.CENTER);
		buttonPanel = new JPanel();
		buttonPanel.setLayout(new GridBagLayout());
		buttonPanel.add(addButton);
		buttonPanel.add(removeButton);
		buttonPanel.add(filterButton);
		JPanel inputsPanel = getInputsPanel();
		JPanel allInputsPanel = new JPanel();
		allInputsPanel.setLayout(new BorderLayout());
		allInputsPanel.add(inputsPanel, BorderLayout.CENTER);
		allInputsPanel.add(buttonPanel, BorderLayout.SOUTH);
		actPanel.add(allInputsPanel, BorderLayout.SOUTH);
		add(selPanel);
		add(actPanel);

	}

	public void init(ObjectTableModel<X> filteredModel, ObjectTableModel<X> selectedModel) {
		init(new JTable(filteredModel), new JTable(selectedModel), filteredModel, selectedModel);
	}

	public abstract JPanel getInputsPanel();

	public abstract void filter();

	public Collection<X> getObjects() {
		ArrayList<X> out = new ArrayList<X>();
		for (int i = 0; i < selectedModel.getRowCount(); i++) {
			out.add(selectedModel.getObject(i));
		}
		return out;
	}

	/*
	 * returns the elements of the filtered list that should be added to the
	 * selected list when the "add" button is pressed
	 */
	public Collection<X> getFilteredForSelected() {
		ArrayList<X> output = new ArrayList<X>();
		int[] inds = filteredList.getSelectedRows();
		for (int i = 0; i < inds.length; i++) {
			X o = filteredModel.getObject(inds[i]);
			output.add(o);
		}
		return output;
	}

	public void add() {
		for (X o : getFilteredForSelected()) {
			if (!selectedModel.contains(o)) {
				selectedModel.addObject(o);
			}
		}
	}

	public void remove() {
		int[] inds = selectedList.getSelectedRows();
		for (int i = inds.length - 1; i >= 0; i--) {
			selectedModel.deleteObject(inds[i]);
		}
	}

	public void filterToSel(MouseEvent e) {
		if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 2) {
			int row = filteredList.rowAtPoint(e.getPoint());
			X x = filteredModel.getObject(row);
			selectedModel.addObject(x);
		}
	}

	public void delSel(MouseEvent e) {
		if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 2) {
			int row = selectedList.rowAtPoint(e.getPoint());
			selectedModel.deleteObject(row);
		}
	}

	public Collection<X> getSelected() {
		return selectedModel.getObjects();
	}

	public int getNumSelected() {
		return selectedModel.getSize();
	}

	public X removeFirstSelected() {
		X val = selectedModel.getObject(0);
		selectedModel.deleteObject(0);
		return val;
	}

	public void addToSelected(Collection<X> objects) {
		if (handlingNewGenome || dataReady) {
			synchronized (asyncSelected) {
				asyncSelected.addAll(objects);
			}
		} else {
			for (X x : objects) {
				selectedModel.addObject(x);
			}
		}

	}

	public void addToSelected(X object) {
		if (handlingNewGenome || dataReady) {
			synchronized (asyncSelected) {
				asyncSelected.add(object);
			}
		} else {
			selectedModel.addObject(object);
		}
	}

	public void close() {
	}

	public boolean isClosed() {
		return true;
	}

}
