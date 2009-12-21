package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;
import javax.swing.JPopupMenu;
import javax.swing.JMenu;
import javax.swing.JMenuItem;

import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.utils.EventSource;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.warpdrive.RectangleLookup;
import edu.mit.csail.cgs.warpdrive.components.*;
import edu.mit.csail.cgs.warpdrive.model.WarpModel;
import edu.mit.csail.cgs.warpdrive.model.ModelProperties;

import java.sql.SQLException;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/* abstract superclass for all Paintables in Warpdrive. 

The WarpDrive framework uses two key classes of objects: WarpPaintables and WarpModels.  The 
model is the datasource.  It queries some underlying data stream (eg, in response to setRegion calls) 
and when the data has been retrieved and processed, it notifies the paintables that use the model.  When
all of a paintable's models are ready, the paintable notifies the enclosing container, which will
call the paint method on the paintable(s). 

Most models and paintables inherit from RegionPaintable and RegionModel to display data across a genomic region.
However, the Model/Paintable system is not specific to that case.

Each model is typically run in a separate thread to improve the overall response time and to prevent the UI
thread from hanging on reltaively long-running database queries.  */

public abstract class WarpPaintable implements VizPaintable, Listener<EventObject>, EventSource<EventObject>, MouseListener {

    /* Listeners are objects that are waiting on the paintable.  When the
       paintable is ready to paint (usually because its underlying model
       has notified the paintable that the underlying model is ready),
       the paintable will notify the listeners that it is ready */
    private HashSet<Listener<EventObject>> listeners;
    private boolean canPaint, wantsPaint;
    /* A common feature of painters is associating String labels with various regions of the display.
       This is the infrastructure for that capability */
    private RectangleLookup<String> labels;
    private int optionKey;
    private Object optionInfo;

    public WarpPaintable () {
        listeners = new HashSet<Listener<EventObject>>();
        labels = new RectangleLookup<String>();
    }

    /* default implementation of Listener<EventObject>.  Notify our listeners
       that we are ready to paint.  */
    public synchronized void notifyListeners(EventObject obj) {
        for (Listener<EventObject> l : listeners) {
            l.eventRegistered(obj);
        }
    }
    
    public void notifyListeners() { 
        notifyListeners(new EventObject(this));
    }
    
    public void addEventListener(Listener<EventObject> l) {
        listeners.add(l);
    }
    public void removeEventListener(Listener<EventObject> l) {
        listeners.remove(l);
    }
    public boolean hasListeners() {return listeners.size() > 0;}

    public abstract PaintableProperties getProperties();
    public boolean canPaint() {return canPaint;}
    public boolean wantsPaint() {return wantsPaint;}
    // only a Paintable should call these methods on itself
    protected void setCanPaint(boolean b) {canPaint = b;}
    protected void setWantsPaint(boolean b) {wantsPaint = b;}
    public void eventRegistered(EventObject e) {}
    public void setLabel(String label) {
        getProperties().TrackLabel = label;
    }
    public String getLabel() {return getProperties().TrackLabel;}
    
    public void cleanup() { 
    }
    
    /** 
     * savePropsInDir will save the painter's properties but
     * subclasses need to take care of saving the model's properties
     * if necessary
     */
    public void savePropsInDir(File dir) {
        PaintableProperties p = getProperties();
        File saveto = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + p.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        p.saveToFile(saveto);
    }
    public void loadPropsInDir(File dir) {
        PaintableProperties p = getProperties();
        File loadfrom = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + p.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        try {
            p.loadFromFile(loadfrom);
        } catch (FileNotFoundException e) {
            System.err.println(e.toString());
        }
    }
    public void saveModelPropsInDir(File dir, WarpModel model) {
        ModelProperties mp = model.getProperties();
        File saveto = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + mp.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        mp.saveToFile(saveto);
    }
    public void loadModelPropsInDir(File dir, WarpModel model) {
        ModelProperties mp = model.getProperties();
        File loadfrom = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + mp.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        mp.saveToFile(loadfrom);
    }    
    /* The order of the keys in this List controls the order of the fields
       in the HashtableConfigurationPanel.  If the List is null, then
       the keys are sorted alphabetically.  Keys present in the hashtable
       but not present in this list are added at the end in an arbitrary order */
    public java.util.List<String> configurationKeyOrder() {return null;}

    /* The label system lets a painter associate rectangles in its
       area with labels.  When the user clicks on the rectangle,
       the label will be displayed.  This relies on the mouseClicked method
       being called (ie, the container must pass the mouse event on */
    /* must be called before addLabel() or clearLabels() */
    public void initLabels() { labels = new RectangleLookup<String>(); }
    
    /* removes all stored labels */
    public void clearLabels(int x, int y, int w, int h, String s) {
        Rectangle rect = new Rectangle(x,y,w,h);
        labels.clearRectangle(rect);
    }    
    public void clearLabels() { labels.clear(); }
    
    /* associates a string with a rectangle.  Clicking on the rectangle
       will show the string */
    public void addLabel(int x, int y, int w, int h, String s) {
        Rectangle rect = new Rectangle(x,y,w,h);
        labels.addValue(s, rect);
    }

    /* returns an ArrayList of JMenuItems (returns null if nothing to be displayed) that
       should be displayed in a popup menu for this mouseclick.
       The container that includes this paintable should gather all of the JMenuItems
       into a single popup menu or other display, preferably labelling them
       by the painter's propertyKey or Label.

       This default implementation will return JMenuItems representing the labels that
       are associated with regions containing the click
    */
    public ArrayList<JMenuItem> mouseClickedMenu(MouseEvent e) {
        ArrayList<JMenuItem> items = new ArrayList<JMenuItem>();
        if (labels != null && e.getButton() == MouseEvent.BUTTON1) {
            int xpos = e.getX();
            int ypos = e.getY();
            final WarpPaintable wp = this;
            boolean any = false;
            Collection<String> labelStrings = labels.getAllValues(new Point(xpos, ypos));
            for(String s : labelStrings) { 
                JMenuItem item = new JMenuItem(s);
                items.add(item);
                item.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                            wp.clickedOnItem(e);
                        }
                    });                
            }
        }
        if (items.size() > 0) {
            return items;
        } else {
            return null;
        }
    }

    /* This is a callback from mouseClickedMenu that allows the painter
       to handle a mouseclick on a menu item that it created 
    */
    public void clickedOnItem(ActionEvent e) {
    }

    /* empty default implementation of MouseListener */
    public void mouseClicked(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}

    /* the option information stores the field (eg bindingScans, bayesresults)
       and value from that field that created this WarpPaintable.  We need to keep
       this information so that a container of WarpPaintables can update a WarpOptions
       structure when a paintable is removed. */ 
    public void setOption(int key, Object info) {
        optionKey = key;
        optionInfo = info;
    }
    public int getOptionKey() {return optionKey;}
    public Object getOptionInfo() {return optionInfo;}

}

