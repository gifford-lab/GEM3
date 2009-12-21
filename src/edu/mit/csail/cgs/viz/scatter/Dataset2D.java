package edu.mit.csail.cgs.viz.scatter;

public class Dataset2D {

    private String labelone, labeltwo;
    private float[][] data;

    /* represents a 2D dataset as a 2xn matrix of floats.  Also stores the label for each
       axis.
    */
    public Dataset2D(float[][] data,
                     String one,
                     String two) {
        this.data = data;
        labelone = one;
        labeltwo = two;
        if (data.length != 2) {
            throw new IllegalArgumentException("data must be a 2xn array");
        }
    }
    public String getLabelOne() {return labelone;}
    public String getLabelTwo() {return labeltwo;}
    public float getVal(int i, int j) {return data[i][j];}
    public int getCount(){return data[0].length;}
}