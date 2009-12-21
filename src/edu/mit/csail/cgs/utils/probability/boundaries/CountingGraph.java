/*
 * Created on Feb 7, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.io.*;
import java.util.*;
import java.text.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.numeric.Numerical;

/**
 * @author tdanford
 */
public class CountingGraph {
    
    private Vector<Layer> layers;
    private int currentN, currentK;

    public CountingGraph(int N, int K) {
        currentN = N; 
        currentK = K;
        layers = new Vector<Layer>();
        Layer nextLayer = null;
        for(int i = 0; i <= N; i++) {
            Layer thisLayer = new Layer(K+1, nextLayer);
            layers.insertElementAt(thisLayer, 0);
            nextLayer = thisLayer;
        }
        
        add();
    }
    
    public double getLogCount() {
        return layers.get(0).getLowestNode().getLogCount(); 
    }
    
    public void arrange(int n, int k, int p0, int s0, boolean equality) { 
        extend(n, k);
        zero();
        setAccessible(false);
        
        int bias = (p0 - s0) + (2 * k) - n;
        //System.out.println("bias: " + bias);
        //System.out.println("n/k: " + n + "," + k + " bias: " + bias);
        
        for(int ki = k; ki >= 0; ki--) { 
            int val = ki * 2;

            for(int ni = 1; ni <= n; ni++) { 
                Node node = layers.get(ni).getNode(ki);
                int upper = ni;
                int lower = ni + bias;
                
                if(equality) {
                    if(lower <= val && val <= upper) { 
                        node.setAccessible(true);
                    }
                } else { 
                    if(lower < val && val < upper) { 
                        node.setAccessible(true);
                    }
                }
            }
            
        }
        
        
        Node firstnode = layers.get(0).getNode(0);
        firstnode.setAccessible(true);
                
        Node lastnode = layers.get(n).getNode(k);
        lastnode.setAccessible(true);
        lastnode.one();
        
        add();        

        //printAccessibility(n, k);
        //printNumbers(n, k);
    }
    
    private void printNumbers(int n, int k) { 
        NumberFormat nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(1);
        
        for(int ki = k; ki >= 0; ki--) {
            for(int ni = 0; ni <= n; ni++) {  
                if(layers.get(ni).getNode(ki).isAccessible()) { 
                    System.out.print(nf.format(Math.exp(layers.get(ni).getNode(ki).getLogCount())) + " ");
                } else { 
                    System.out.print("X ");
                }
            }
            System.out.println();
        }        
    }
    
    private void printAccessibility(int n, int k) { 
        for(int ki = k; ki >= 0; ki--) {
            for(int ni = 0; ni <= n; ni++) {  
                if(layers.get(ni).getNode(ki).isAccessible()) { 
                    System.out.print("*");
                } else { 
                    System.out.print("-");
                }
            }
            System.out.println();
        }
    }
    
    public void extend(int N, int K) {
        if(N > currentN || K > currentK) { 
            int extra = N+1-layers.size();
            
            for(int i = layers.size()-1; i >= 0; i--) { 
                layers.get(i).extend(K+1);
            }
            
            Layer nextLayer = layers.get(0);
            for(int i = extra-1; i >= 0; i--) {
                Layer thisLayer = new Layer(K+1, nextLayer);
                layers.add(0, thisLayer);
                nextLayer = thisLayer;
            }
            
            currentN = N; 
            currentK = K;
            
            add();
        }
    }
    
    private void setAccessible(boolean v) { 
        for(Layer l : layers) { l.setAccessible(v); }
    }
    
    private void add() { 
        for(int i = layers.size()-1; i >= 0; i--) { 
            layers.get(i).add();
        }
    }
    
    private void zero() { 
        for(Layer l : layers) { l.zero(); }
    }
    
    public int getN() { return currentN; }
    public int getK() { return currentK; }
    
    public Vector<String> getGridStrings() { 
        Vector<String> strv = new Vector<String>();
        Vector<char[]> layerStrings = new Vector<char[]>();
        Vector<StringBuilder> sbv = new Vector<StringBuilder>();
        
        for(Layer l : layers) { layerStrings.add(l.getGridString()); }
        
        for(int i = 0; i <= currentK; i++) {
            StringBuilder sb = new StringBuilder();
            sbv.add(sb);
            for(char[] array : layerStrings) { 
                sb.append(array[i]);
            }
        }

        for(StringBuilder sb : sbv) { strv.add(sb.toString()); }
        return strv;
    }
    
    private static class Layer { 
        
        private Vector<Node> nodes;
        private Layer nextLayer;
        
        public Layer(int height, Layer next) {
            nodes = new Vector<Node>();
            
            for(int i = 0; i < height; i++) {
                Node newNode = new Node();
                nodes.add(newNode);
            }
            
            setNextLayer(next);
        }
        
        public void setAccessible(boolean v) { 
            for(Node n : nodes) { n.setAccessible(v); }
        }
        
        public void setNextLayer(Layer l) { 
            nextLayer = l;
            for(int i = 0; i < nodes.size(); i++) { 
                setNodeNeighbors(i);
            }
        }
        
        private void setNodeNeighbors(int nodeIndex) {
            Node n = nodes.get(nodeIndex);
            Node lower = null, upper = null;
            if(nextLayer != null) { 
                lower = nextLayer.nodes.get(nodeIndex);
                if(nodeIndex+1 < nextLayer.nodes.size()) { 
                    upper = nextLayer.nodes.get(nodeIndex+1);
                }
                n.setLowerNode(lower);
                n.setUpperNode(upper);
            }
        }
        
        public void extend(int newHeight) { 
            while(nodes.size() < newHeight) { 
                nodes.add(new Node()); 
            }
            
            for(int i = 0; i < nodes.size(); i++ ) {
                setNodeNeighbors(i);
            }
        }
        
        public Node getLowestNode() { return nodes.get(0); }
        public Node getHighestNode() { return nodes.get(nodes.size()-1); }
        
        public Node getNode(int i) {
            if(i < 0 || i >= nodes.size()) { 
                throw new IllegalArgumentException("Can't select #" + i + " from " + nodes.size());
            }
            return nodes.get(i); 
        }
        
        public void zero() { 
            for(int i = 0; i < nodes.size(); i++) { 
                nodes.get(i).zero();
            }
        }
        
        public void add() { 
            for(int i = 0; i < nodes.size(); i++) { 
                nodes.get(i).add();
            }
        }
        
        public char[] getGridString() { 
            char[] array = new char[nodes.size()];
            for(int i = 0; i < array.length; i++) { 
                if(nodes.get(i).isAccessible()) { 
                    array[i] = '*';
                } else { 
                    array[i] = '-';
                }
            }
            return array;
        }
    }
    
    private static class Node {
        
        private boolean accessible;
        private Double logCount;
        private Node nextLower, nextUpper;
        
        public Node() {
            nextLower = nextUpper = null;
            logCount = null;
            accessible = true;
        }
        
        public void setLowerNode(Node l) { nextLower = l; }
        public void setUpperNode(Node u) { nextUpper = u; }
        
        public double getLogCount() {
            if(logCount == null) { return -Double.MAX_VALUE; }
            return logCount; 
        }
        
        public void setAccessible(boolean v) {
            accessible = v;
        }
        
        public boolean isAccessible() { return accessible; }
        
        public void one() { 
            logCount = 0.0;
        }
        
        public void setCount(int c) { 
            logCount = Math.log((double)c);
        }
        
        public void zero() { logCount = null; }
        public boolean isZero() { return logCount == null; }
        
        public void add() {
            if(accessible) { 
                if(nextLower != null && nextUpper != null) {
                    if(!nextLower.isZero() && !nextUpper.isZero()) { 
                        logCount = Numerical.log_add(nextLower.logCount, nextUpper.logCount);
                    } else { 
                        if(nextLower.isZero()) { logCount = nextUpper.logCount; }
                        if(nextUpper.isZero()) { logCount = nextLower.logCount; }
                    }
                } else { 
                    if(nextLower != null) { logCount = nextLower.logCount; }
                    if(nextUpper != null) { logCount = nextUpper.logCount; }
                }
            }
        }
    }
    
    private void printStrings(Collection<String> strs) { printStrings(strs, System.out); } 
    private void printStrings(Collection<String> strs, PrintStream ps) { for(String s : strs) { ps.println(s); } }
}
