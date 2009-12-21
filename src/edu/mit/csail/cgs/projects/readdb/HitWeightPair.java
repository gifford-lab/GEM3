package edu.mit.csail.cgs.projects.readdb;

public class HitWeightPair implements Comparable<HitWeightPair> {
    public int pos;
    public float weight;
    public HitWeightPair (int p, float w) {
        pos = p;
        weight = w;
    }
    public boolean equals(Object o) {
        if (o instanceof HitWeightPair) {
            HitWeightPair other = (HitWeightPair)o;
            return pos == other.pos && weight == other.weight;
        } else {
            return false;
        }
    }
    public int compareTo(HitWeightPair o) {
        return pos - o.pos;
    }
}