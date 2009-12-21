package edu.mit.csail.cgs.datasets.chipchip;

public class ExptNameVersionChannel extends ExptNameVersion {

    public static final int IP = 1, WCE = 2;

    private int channel;

    public ExptNameVersionChannel(String n, String v, String r, int c) {
        super(n,v,r);
        if (c != IP && c != WCE) {
            throw new IllegalArgumentException("Unknown channel : " + c);
        }
        channel = c;
    }
    public ExptNameVersionChannel(String n, String v, int c) {
        super(n,v);
        if (c != IP && c != WCE) {
            throw new IllegalArgumentException("Unknown channel : " + c);
        }
        channel = c;
    }
    public ExptNameVersionChannel(String pieces[]) {
        super(pieces[1],pieces[1]);
        if (pieces.length == 3) {
            channel = Integer.parseInt(pieces[2]);
            if (channel != IP && channel != WCE) {
                throw new IllegalArgumentException("Unknown channel : " + channel);
            }   
        } else if (pieces.length == 4) {
            setReplicate(pieces[2]);
            channel = Integer.parseInt(pieces[3]);
            if (channel != IP && channel != WCE) {
                throw new IllegalArgumentException("Unknown channel : " + channel);
            }   
        } else {
            throw new IllegalArgumentException("Can't process " + pieces + " into an ENVC");
        }



    }

    public int getChannel() {return channel;}

}
