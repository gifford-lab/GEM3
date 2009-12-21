package edu.mit.csail.cgs.ewok.verbs;

public class Counter<X extends Object> implements Distiller<X,Integer> {
    private int count;
    public Counter() {
        count = 0;
    }
    public Integer execute(X o) {
        return count++;
    }
    public Integer getCurrent() {
        return count;
    }
    public void reset() {
        count = 0;
    }
}
