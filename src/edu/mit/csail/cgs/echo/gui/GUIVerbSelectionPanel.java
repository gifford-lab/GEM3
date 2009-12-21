/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

import edu.mit.csail.cgs.echo.*;
import edu.mit.csail.cgs.utils.EventSource;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.viz.components.SelectionEvent;

import edu.mit.csail.cgs.echo.components.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;

/**
 * @author tdanford
 */
public class GUIVerbSelectionPanel extends JPanel implements EventSource<SelectionEvent<GUIVerbSelector>> {
    
    public static void main(String[] args) {
        System.out.println("Startup.");
        
        JFrame f = new JFrame("Test Window");
        Container c = f.getContentPane();
        c.setLayout(new BorderLayout());
        GUIVerbSelectionPanel p = new GUIVerbSelectionPanel();
        c.add(p, BorderLayout.CENTER);
        
        p.addVerb(new GUIVerbSelector(new ListingSink()));
        p.addVerb(new GUIVerbSelector(new ChromRegionWrapper()));
        p.addVerb(new GUIVerbSelector(new RefGeneGenerator()));
        
        f.setVisible(true);
        f.pack();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        p.addEventListener(new Listener<SelectionEvent<GUIVerbSelector>> () {
            public void eventRegistered(SelectionEvent<GUIVerbSelector> e) {
                System.out.println("Selected: " + e.getFirstValue().getName());
            } 
        });
    }
    
    private Vector<GUIVerbSelector> selectors;
    private Set<ClassSelector> inputClassSelectors, outputClassSelectors;
    
    private DefaultComboBoxModel inputModel, outputModel;
    private JComboBox inputBox, outputBox;
    
    private DefaultListModel selectorModel;
    private JList selectorList;
    
    private JButton filterButton;
    private ClassSelector none, any;
    
    private JTextField paramField;
    
    private EventSource.Default<SelectionEvent<GUIVerbSelector>> src;
    
    public GUIVerbSelectionPanel() { 
        super();
        
        src = new EventSource.Default<SelectionEvent<GUIVerbSelector>>();
        
        none = new ClassSelector(NO_CLASS);
        any = new ClassSelector(ANY_CLASS);
        
        selectors = new Vector<GUIVerbSelector>();
        inputClassSelectors = new HashSet<ClassSelector>();
        outputClassSelectors = new HashSet<ClassSelector>();
        
        inputModel = new DefaultComboBoxModel();
        outputModel = new DefaultComboBoxModel();
        
        inputModel.addElement(none); inputModel.addElement(any);
        outputModel.addElement(none); outputModel.addElement(any);

        inputBox = new JComboBox(inputModel);
        outputBox = new JComboBox(outputModel);
        
        paramField = new JTextField();
        
        inputModel.setSelectedItem(any);
        outputModel.setSelectedItem(any);
        
        selectorModel = new DefaultListModel();
        selectorList = new JList(selectorModel);

        setLayout(new BorderLayout());

        JPanel fromPanel = new JPanel(); 
        fromPanel.setLayout(new BorderLayout());
        fromPanel.setBorder(new TitledBorder("From:"));
        fromPanel.add(inputBox, BorderLayout.NORTH);
        
        JPanel toPanel = new JPanel();
        toPanel.setLayout(new BorderLayout());
        toPanel.setBorder(new TitledBorder("To:"));
        toPanel.add(outputBox, BorderLayout.NORTH);
        
        JPanel topPanel = new JPanel();
        topPanel.setLayout(new GridLayout(1, 2));
        topPanel.add(fromPanel); 
        topPanel.add(toPanel);
        add(topPanel, BorderLayout.NORTH);
        
        JPanel listPanel = new JPanel();
        listPanel.setLayout(new BorderLayout());
        listPanel.setBorder(new TitledBorder("Available:"));
        listPanel.add(new JScrollPane(selectorList));
        add(listPanel, BorderLayout.CENTER);
        
        JPanel bottomPanel = new JPanel();
        bottomPanel.setLayout(new GridLayout(1, 1));
        bottomPanel.add(paramField);
        bottomPanel.setBorder(new TitledBorder("Parameters:"));
        add(bottomPanel, BorderLayout.SOUTH);
        
        inputBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) { 
                filter();
            }
        });
        
        outputBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) { 
                filter();
            }
        });
        
        selectorList.addMouseListener(new MouseAdapter() { 
            public void mouseClicked(MouseEvent e) { 
                if(e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 2) { 
                    int index = selectorList.locationToIndex(e.getPoint());
                    GUIVerbSelector sel = (GUIVerbSelector)selectorModel.elementAt(index);
                    sel = new GUIVerbSelector(sel, paramField.getText());
                    SelectionEvent<GUIVerbSelector> evt = 
                        new SelectionEvent<GUIVerbSelector>(GUIVerbSelectionPanel.this, sel);
                    src.fireEvent(evt);
               }
            }
        });
        
        filter();
    }
    
    public void addEventListener(Listener<SelectionEvent<GUIVerbSelector>> el) {
        src.addEventListener(el);
    }

    public boolean hasListeners() {
        return src.hasListeners();
    }

    public void removeEventListener(Listener<SelectionEvent<GUIVerbSelector>> el) {
        src.removeEventListener(el);
    }

    public void addVerb(GUIVerbSelector s) { 
        if(!selectors.contains(s)) { 
            selectors.add(s);
            
            SelfDescribingVerb v = s.getVerb();
            EchoType[] inputClasses = v.getInputClasses();
            for(int i = 0; inputClasses != null && i < inputClasses.length; i++) { 
                ClassSelector cs = new ClassSelector(inputClasses[i]);
                if(!inputClassSelectors.contains(cs)) { 
                    inputClassSelectors.add(cs);
                    inputModel.addElement(cs);
                }
            }
            
            EchoType outputClass = v.getOutputClass();
            if(outputClass != null) { 
                ClassSelector cs = new ClassSelector(outputClass);
                if(!outputClassSelectors.contains(cs)) {
                    outputClassSelectors.add(cs);
                    outputModel.addElement(cs);
                }
            }
            
            filter();
        }
    }

    public void filter() { 
        ClassSelector fromsel = (ClassSelector)inputModel.getSelectedItem();
        ClassSelector tosel = (ClassSelector)outputModel.getSelectedItem();

        selectorModel.clear();
        for(GUIVerbSelector sel : selectors) { 
            if(matchesSelectors(sel, fromsel, tosel)) { 
                selectorModel.addElement(sel);
            }
        }
    }
    
    private static boolean matchesSelectors(GUIVerbSelector verbsel, ClassSelector insel, ClassSelector outsel) { 
        SelfDescribingVerb v = verbsel.getVerb();

        if(insel.type != ANY_CLASS) { 
            EchoType[] inputClasses = v.getInputClasses();
            if(inputClasses == null) { 
                if(insel.type == SPECIFIC_CLASS) { return false; }
            } else { 
                if(insel.type == NO_CLASS) { return false; }
                boolean foundMatch = false;
                for(int i = 0; !foundMatch && i < inputClasses.length; i++) { 
                	foundMatch = insel.selected.isSubType(inputClasses[i]);
                }
                if(!foundMatch) { return false; }
            }
        }
        

        if(outsel.type != ANY_CLASS) { 
            EchoType outputClass = v.getOutputClass();
            if(outputClass == null) { 
                if(outsel.type == SPECIFIC_CLASS) { return false; }
            } else { 
                if(outsel.type == NO_CLASS) { return false; }
                if(outputClass.isSubType(outsel.selected)) { return false; }
            }
        }
        
        return true;
    }
    
    private static final int ANY_CLASS = 0;
    private static final int NO_CLASS = 1;
    private static final int SPECIFIC_CLASS = 2;
    
    private static class ClassSelector {
        
        public int type;
        public String name;
        public EchoType selected;
        
        public ClassSelector(int t) { 
            type = t;
            if(type != ANY_CLASS && type != NO_CLASS) { 
                throw new IllegalArgumentException(); 
            }
            
            if(type == ANY_CLASS) { name = "ANY"; }
            if(type == NO_CLASS) { name = "NONE"; }
            selected = null;
        }
        
        public ClassSelector(EchoType s) { 
            type = SPECIFIC_CLASS;
            selected = s;
            name = s.getName();
            
            /*
            String cn = s.getName();
            String[] array = cn.split("\\.");
            name = array.length > 0 ? array[array.length-1] : cn;
            */
        }
        
        public String toString() { return name; }
        
        public boolean equals(Object o) { 
            if(!(o instanceof ClassSelector)) { return false; }
            ClassSelector cs = (ClassSelector)o;
            if(type != cs.type) { return false; }
            if(!name.equals(cs.name)) { return false; }
            if(cs.selected != null && selected == null) { return false; }
            if(selected != null && cs.selected == null) { return false; }
            if(selected != null && cs.selected != null) { 
                if(!selected.equals(cs.selected)) { return false; }
            }
            return true;
        }
        
        public int hashCode() { 
            int code = 17;
            code += type; code *= 37;
            code += name.hashCode(); code *= 37;
            return code;
        }
    }

}

