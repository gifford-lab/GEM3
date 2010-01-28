package edu.mit.csail.cgs.projects.readdb;

import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;

class Request {
    public String type, alignid;
    public Integer chromid, start, end;
    public Boolean isPaired, isLeft, isPlusStrand;
    public Float minWeight;
    public Map<String,String> map;
    public List<String> list;
    public Request () {
        map = new HashMap<String,String>();
        list = new ArrayList<String>();
        type = null;
        alignid = null;
        chromid = null;
        start = null;
        end = null;
        isPaired = null;
        isLeft = null;
        isPlusStrand = null;
        minWeight = null;
    }
    public void clear() {
        type = null;
        alignid = null;
        chromid = null;
        start = null;
        end = null;
        isPaired = null;
        isLeft = null;
        isPlusStrand = null;
        minWeight = null;
        map.clear();
        list.clear();
    }
    /**
     * parses from a list of strings in the same format as toString() outputs.
     * returns null on success or an error message on failure.
     */
    public String parse(List<String> args) {
        clear();
        for (String s : args) {
            String pieces[] = s.split("\\s*=\\s*");
            if (pieces.length == 2) {
                if (pieces[0].equals("alignid")) {
                    alignid = pieces[1];
                } else if (pieces[0].equals("requesttype")) {
                    type = pieces[1];
                } else if (pieces[0].equals("chromid")) {
                    try {
                        chromid = new Integer(pieces[1]);
                    } catch (NumberFormatException e) {
                        return "invalid number for chromid " + pieces[1] ;
                    }
                } else if (pieces[0].equals("start")) {
                    try {
                        start = new Integer(pieces[1]);
                    } catch (NumberFormatException e) {
                        return "invalid number for start " + pieces[1] ;
                    }
                } else if (pieces[0].equals("end")) {
                    try {
                        end = new Integer(pieces[1]);
                    } catch (NumberFormatException e) {
                        return "invalid number for end " + pieces[1] ;
                    }
                } else if (pieces[0].equals("ispaired")) {
                    isPaired = new Boolean(pieces[1]);
                } else if (pieces[0].equals("isleft")) {
                    isLeft = new Boolean(pieces[1]);
                } else if (pieces[0].equals("isplusstrand")) {
                    isPlusStrand = new Boolean(pieces[1]);
                } else if (pieces[0].equals("minweight")) {
                    try {
                        minWeight = new Float(pieces[1]);
                    } catch (NumberFormatException e) {
                        return "invalid minweight " + pieces[1] ;
                    }
                } else {
                    map.put(pieces[0],pieces[1]);
                }
            } else if (pieces.length == 1) {
                list.add(s);
            } else {
                return "Invalid number of fields on line " + s;
            }
        }
        if (type == null) {
            return "must provide a requestype";
        }
        if (isPaired && isLeft == null) {
            return "must provide isleft when providing ispaired";
        }

        return null;
    }
    public String toString() {
        StringBuffer out = new StringBuffer();
        if (type != null) {
            out.append("requesttype=" + type + "\n");
        }
        if (alignid != null) {
            out.append("alignid=" + alignid + "\n");
        }
        if (chromid != null) {
            out.append("chromid=" + chromid + "\n");
        }
        if (start != null) {
            out.append("start=" + start + "\n");            
        }
        if (end != null) {
            out.append("end=" + end + "\n");
        }
        if (isPaired != null) {
            out.append("ispaired=" + isPaired + "\n");
        }
        if (isLeft != null) {
            out.append("isleft=" + isLeft + "\n");
        }
        if (isPlusStrand != null) {
            out.append("isplusstrand=" + isPlusStrand + "\n");
        }
        for (String k : map.keySet()) {
            out.append(k + "=" + map.get(k) + "\n");
        }
        for (String l : list) {
            out.append(l + "\n");
        }
        out.append("ENDREQUEST\n");
        return out.toString();
    }
}