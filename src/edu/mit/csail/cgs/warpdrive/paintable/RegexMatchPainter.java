package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionMapperModel;

public class RegexMatchPainter extends RegionPaintable {
    
    private RegionMapperModel<String> model;
    private ColorSet cs;
    private RegexMatchProperties props;

    public RegexMatchPainter(RegionMapperModel<String> model) {
        super();
        this.model = model;
        model.addEventListener(this);
        cs = new ColorSet();
        props = new RegexMatchProperties(10);
    }
    public RegexMatchProperties getProperties() {return props;}

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }
    public void addRegex(String r) {
        props.addRegex(r);
    }
    public void addRegex(String label, String regex) {
        props.addRegex(label,regex);
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public synchronized void eventRegistered(EventObject e) {
        if (e.getSource() == model &&
            model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }    

    public int getMaxVertSpace() { 
        return 40;
    }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        if (!canPaint()) {
            return;
        }
        int w = x2 - x1;
        int h = y2 - y1;
        Region region = getRegion();
        int regionstart = region.getStart(), regionend = region.getEnd();
        int regionwidth = regionend - regionstart;
        String wholestring = model.getResults();
        char[] chararray = wholestring.toCharArray();
        int length = chararray.length;
        SequenceUtils.reverseComplement(chararray);
        String revcomp = new String(chararray);

        if (wholestring == null) {
            return;
        }
        cs.reset();
        clearLabels();
        ArrayList<RegexHit> hits = new ArrayList<RegexHit>();
        
        for (int i = 0; i < props.expressions.length; i++) {
            if (props.expressions[i] == null) {
                continue;
            }
            if (props.expressions[i].length() == 0) {
                continue;
            }

            String label = props.labels[i];
            Pattern pattern = Pattern.compile(props.expressions[i]);
            if (pattern == null) {continue;}
            Matcher matcher = pattern.matcher(wholestring);
            while (matcher.find()) {
                String chars = matcher.group();
                int start = matcher.start();
                int end = matcher.end();
                RegexHit hit = new RegexHit(start,end,label,true,false,cs.getColor(chars));
                hit.labels.add(chars);
                if (hits.contains(hit)) {
                    hits.get(hits.indexOf(hit)).combine(hit);
                } else {
                    hits.add(hit);
                }
            }
            matcher = pattern.matcher(revcomp);
            while (matcher.find()) {
                String chars = matcher.group();
                int end = length - matcher.start();
                int start = length - matcher.end();
                RegexHit hit = new RegexHit(start,end,label,false,true,cs.getColor(chars));
                hit.labels.add(chars);
                if (hits.contains(hit)) {
                    hits.get(hits.indexOf(hit)).combine(hit);
                } else {
                    hits.add(hit);
                }
            }
        }

        for (RegexHit hit : hits) {
            int startx = getXPos(hit.start, 0, regionwidth, x1, x2);
            int endx = getXPos(hit.end, 0, regionwidth, x1, x2);
            g.setColor(hit.color);
            g.drawRect(startx, y1+1, Math.max(endx - startx,1), h-2);
            if (hit.plus) {
                g.drawLine(startx,y1+1, endx, y1 + h/2);
                g.drawLine(startx,y1+h, endx, y1 + h/2);
            }
            if (hit.minus) {
                g.drawLine(startx, y1 + h/2, endx, y1+1);
                g.drawLine(startx, y1 + h/2, endx, y1+h);
            }

            for (String l : hit.labels) {
                addLabel(startx, y1, endx - startx, h, l);            
            }
        }
    }     
}

class RegexHit {

    boolean plus, minus;
    public int start, end;
    public Set<String> labels;
    public Color color;
    public RegexHit (int start, int end, String label, boolean plus, boolean minus, Color color) {
        this.labels = new HashSet<String>();        
        this.labels.add(label);
        this.start = start;
        this.end = end;
        this.plus = plus;
        this.minus = minus;
        this.color = color;
    }
    public int Hashcode() {
        return start * 5 + end * 6;
    }
    public void combine(RegexHit other) {
        plus = plus || other.plus;
        minus = minus || other.minus;
        labels.addAll(other.labels);
    }
    public boolean equals(Object o) {
        if (o instanceof RegexHit) {
            RegexHit r = (RegexHit) o;
            return r.start == start && r.end == end;
        }
        return false;
    }
}

