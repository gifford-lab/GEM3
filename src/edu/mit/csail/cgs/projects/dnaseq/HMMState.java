package edu.mit.csail.cgs.projects.dnaseq;

class HMMState {
    
    private int A, C, T, G;
    private int[] readCounts;
    private long totalObs;
    public HMMState() {
        A = 0;
        C = 0;
        T = 0;
        G = 0;
        readCounts = new int[10];
    }
    public HMMState(int a, int c, int g, int t, int[] counts) {
        A = a;
        C = c;
        T = t;
        G = g;
        readCounts = counts;
        totalObs = A + C + T + G;
    }
    public void addData(char letter,
                        int reads) {
        if (letter == 'A' || letter == 'a') {
            A++;
        } else if (letter == 'C' || letter == 'c') {
            C++;
        } else if (letter == 'G' || letter == 'g') {
            G++;
        } else if (letter == 'T' || letter == 't') {
            T++;
        }
        if (reads >= readCounts.length) {
            reads = readCounts.length - 1;
        }
        readCounts[reads]++;
        totalObs++;
    }
    public int[] getCounts() {return readCounts;}
    public double getProb(char letter,
                          int reads) {

        double lcount, rcount;
        if (letter == 'A' || letter == 'a') {
            lcount = A;
        } else if (letter == 'C' || letter == 'c') {
            lcount = C;
        } else if (letter == 'G' || letter == 'g') {
            lcount = G;
        } else if (letter == 'T' || letter == 't') {
            lcount = T;
        } else {
            lcount = 1;
        }

        if (reads >= readCounts.length) {
            reads = readCounts.length - 1;
        }
        rcount = readCounts[reads];
        
        return (lcount / totalObs) * (rcount / totalObs);

    }
    public String toString() {
        double avg = 0;
        for (int i = 1; i < readCounts.length; i++) {
            avg += readCounts[i] * i;
        }
        avg = avg / totalObs;

        StringBuilder sb = new StringBuilder(String.format("A=%d C=%d G=%d T=%d  avg=%.2f",A,C,G,T,avg));
        for (int i = 0; i < readCounts.length; i++) {
            sb.append("\t" +readCounts[i]);
        }
        return sb.toString();
    }
    public String serialize() {
        StringBuilder sb = new StringBuilder(String.format("%d\t%d\t%d\t%d",A,C,G,T));
        for (int i = 0; i < readCounts.length; i++) {
            sb.append("\t"+readCounts[i]);
        }
        return sb.toString();
    }
    public static HMMState deserialize(String string) {
        HMMState s = new HMMState();
        String pieces[] = string.split("\\t");
        s.A = Integer.parseInt(pieces[0]);
        s.C = Integer.parseInt(pieces[1]);
        s.G = Integer.parseInt(pieces[2]);
        s.T = Integer.parseInt(pieces[3]);
        s.totalObs = s.A + s.C + s.G + s.T;
        for (int i = 0; i < s.readCounts.length; i++) {
            s.readCounts[i] = Integer.parseInt(pieces[4+i]);
        }
        return s;
    }
}