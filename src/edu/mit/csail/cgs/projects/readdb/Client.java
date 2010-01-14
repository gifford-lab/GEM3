package edu.mit.csail.cgs.projects.readdb;

import java.net.*;
import java.util.*;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import javax.security.sasl.*;
import javax.security.auth.callback.*;

/**
 * API for remote access to the readdb server.
 *
 *  Calls throw IOException on network errors.
 *  Calls throw ClientException on other errors (authentication, authorization, invalid request ,etc.
 *
 *  Note that the Client doesn't know anything about strandedness.  If you want to distinguish
 *  between hits on the + and - strands, you need to have one chromosome for each strand, eg, 1+ and 1-.  
 *  This is the default format that ImportHits uses.
 *
 *  Client generally assumes, but doesn't require, that the hit positions are the 5' end of the
 *  hit.  You can do something else, though it'll break the extended histogramming
 *
 * Client IS NOT REENTRANT.  Do not overlap calls to a single Client object.
 * 
 */
public class Client implements ReadOnlyClient {    

    public static final String SaslMechanisms[] = {"CRAM-MD5","DIGEST-MD5"};
    /* socket is the socket to talk to the server.  outstream and instream are from
       the socket
        */       
    Socket socket;
    OutputStream outstream;
    BufferedInputStream instream;
    /* temporary space for receiving data; contents not persistent between method calls */
    byte[] buffer;
    private static final int BUFFERLEN = 8192*20;
    ByteOrder myorder;


    public Client (String hostname,
                   int portnum,
                   String username,
                   String passwd) throws IOException, ClientException {
        init(hostname,portnum,username,passwd);
    }
    /**
     * creates the default connection
     * as specified by ~/.readdb_passwd or a readdb_passwd found in the classpath
     *
     * Must have keys hostname, port, username, and passwd in a format
     * that java.util.PropertyResourceBundle can read
     */
    public Client() throws IOException, ClientException {
        String homedir = (String)System.getProperties().get("user.home");
        String basename = "readdb_passwd";
        if (System.getenv("READDBROLE") != null) {
            basename = System.getenv("READDBROLE") + basename;
        }
        String fname = homedir + "/." + basename;
        File propfile = new File(fname);
        PropertyResourceBundle bundle = null;
        if (propfile.exists() && propfile.canRead()) {
            bundle = new PropertyResourceBundle(new FileInputStream(propfile));
        } else {
            ClassLoader cl = ClassLoader.getSystemClassLoader();
            URL url = cl.getResource(basename);
            if (url != null) {
                bundle = new PropertyResourceBundle(url.openStream());
            } else {
                throw new IOException("Can't read connection properties from " + fname);
            }
        }
        String hostname = bundle.getString("hostname");
        String port = bundle.getString("port");
        String username = bundle.getString("username");
        String password = bundle.getString("passwd");
        init(hostname, Integer.parseInt(port), username, password);
    }
    
    private void init(String hostname,
                   int portnum,
                   String username,
                   String passwd) throws IOException, ClientException {
        myorder = ByteOrder.nativeOrder();
        socket = new Socket(hostname,portnum);
        socket.setTcpNoDelay(true);
        socket.setSendBufferSize(BUFFERLEN);
        socket.setReceiveBufferSize(BUFFERLEN);
        /* linger = true, time = 0 means that the other side gets a reset if the
           socket is closed on our end (eg, java exits).  We turn this off just before
           sending a "bye" to allow for a graceful shutdown.  But the RST in
           other cases lets the server figure out that we've disappeared
        */
        socket.setSoLinger(true,0);
        outstream = socket.getOutputStream();
        instream = new BufferedInputStream(socket.getInputStream());
        buffer = new byte[BUFFERLEN];
        if (!authenticate(hostname,username,passwd)) {
            throw new ClientException("Authentication Exception Failed");
        }
        sendString("byteorder\n" + (ByteOrder.nativeOrder() == ByteOrder.BIG_ENDIAN ? "big" : "little") + "\nENDREQUEST\n");
        String output = readLine();
        if (!output.equals("OK")) {
            throw new ClientException("Couldn't set byte order " + output);
        }
    }
    
    /**
     * performs the SASL authentication exchange with the server.  currently called by the constructor
     */
    private boolean authenticate(String hostname, String username, String password) throws IOException {
        SaslClient sasl = null;
        try {
            sendString(username + "\n");
            Map<String,String> props = new HashMap<String,String>();
            props.put("Sasl.POLICY_NOPLAINTEXT","true");
            props.put("Sasl.POLICY_NOANONYMOUS","true");
            sasl = Sasl.createSaslClient(SaslMechanisms,
                                         username,
                                         "readdb",
                                         hostname,
                                         props,
                                         new ClientCallbackHandler(username,password));
            if (sasl == null) { return false; }
            byte[] response = (sasl.hasInitialResponse() ? sasl.evaluateChallenge(new byte[0]) : new byte[0]);
            byte continueByte = 1;
            while (continueByte != 0) {
                outstream.write((response.length + "\n").getBytes());
                outstream.write(response);
                outstream.flush();
                //                System.err.println("Sent " + response.length + " bytes");
                int length = Integer.parseInt(readLine());
                //                System.err.println("Waiting to read " + length + " bytes");
                byte[] challenge = new byte[length];
                int read = 0;
                while (read < length) {
                    read += instream.read(challenge, read, length - read);
                    //                    System.err.println("   read " + read);
                }
                /* the continueByte tells us whether the server
                   expects to do another round.  Necessary because sometimes
                   isComplete() returned true here but the server wasn't done
                */
                continueByte = (byte)instream.read();
                //                System.err.println("Continue byte is " + continueByte);
                if (!sasl.isComplete()) { 
                    response = sasl.evaluateChallenge(challenge);
                } else {
                    response = new byte[0];
                }
            }
            sasl.dispose();
            String status = readLine();
            return (status.equals("authenticated as " + username));
        } catch (SaslException e) {
            e.printStackTrace();
            if (sasl != null) {
                sasl.dispose();
            }
            return false;
        }
    }
    
    /** sends a string to the server and flushes the socket 
     */
    private void sendString(String s) throws IOException {
        //        System.err.println("SENDING " + s);
        outstream.write(s.getBytes());
        outstream.flush();
    }
    /** sends a string to the server
     */
    private void sendStringNoFlush(String s) throws IOException {
        outstream.write(s.getBytes());
    }
    /** reads one line from the server.  blocking.
     */
    private String readLine() throws IOException {
        int pos = 0;
        int i;
        while ((i = instream.read()) != -1) {
            if (i == '\n') {
                break;
            } else {
                buffer[pos++] = (byte)i;
            }
        }
        String out = new String(buffer,0,pos);
        //        System.err.println("READ " + out);
        return out;
    }
    public void shutdown() throws IOException, ClientException{
        sendString("shutdown\nENDREQUEST\n");
    }
    /** Stores hits for the specified alignment and chromosome.
     *  Will append to the current set of hits if they exist.
     *  Throws ClientException if the currently authenticated user doesn't
     *  have permission to create the alignment when it doesn't already exist.
     */
    public void store(String alignid, String chromid, int hits[], float weights[]) throws IOException, ClientException {
        sendString("store\n" + alignid + "\n" + chromid + "\n" + hits.length + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        Bits.sendInts(hits, outstream,buffer, myorder);
        Bits.sendFloats(weights, outstream,buffer, myorder);
        System.err.println("Sent " + hits.length + " hits to the server");
        outstream.flush();        
        response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
    }
    /** Returns true if the alignment and chromosome exist and are accessible
     * to the user.  Returns false if they don't exist or if they aren't accessible
     */
    public boolean exists(String alignid) throws IOException {
        sendString("exists\n" + alignid + "\nENDREQUEST\n");
        String response = readLine();
        if (response.equals("exists")) {
            return true;
        } else if (response.equals("unknown")) {
            return false;
        } else {
            return false;
        }
    }
    /**
     * Deletes an alignment (all chromosomes and the acl).
     */
    public void deleteAlignment(String alignid) throws IOException, ClientException {
        sendString("deletealign\n" + alignid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
    }
    /**
     * Returns the set of chromosomes that exist for this alignment
     */
    public Set<String> getChroms(String alignid) throws IOException, ClientException {
        sendString("getchroms\n" + alignid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numchroms = Integer.parseInt(readLine());
        Set<String> output = new HashSet<String>();
        while (numchroms-- > 0) {
            output.add(readLine());
        }
        return output;
    }
    /**
     * Returns the total number of hits in this alignment
     */
    public int getCount(String alignid) throws IOException, ClientException {
        int count = 0;
        for (String c : getChroms(alignid)) {
            count += getCount(alignid, c);
        }
        return count;
    }
    /**
     * Returns the sum of the weights of all hits in this alignment
     */
    public double getWeight(String alignid) throws IOException, ClientException {
        double total = 0;
        for (String c : getChroms(alignid)) {
            total += getWeight(alignid, c);
        }
        return total;
    }

    /** returns the total number of hits on the specified chromosome in the alignment
     */
    public int getCount(String alignid, String chromid)  throws IOException, ClientException {
        sendString("count" + "\n" + alignid + "\n" + chromid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return numhits;
    }
    /** returns the total weight on the specified chromosome in this alignment
     */
    public double getWeight(String alignid, String chromid) throws IOException, ClientException {
        sendString("weight" + "\n" + alignid + "\n" + chromid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        return Double.parseDouble(readLine());
    }
    /** returns the number of hits in some range of chromosomal coordinates in the alignment
     */
    public int getCountRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException {
        return getCountRange(alignid, chromid, start, stop, Float.NaN);
    }
    /** returns the number of hits in some range of chromosomal coordinates in the alignment.  Only
     *  includes hits whose weight is greater than minweight
     */
    public int getCountRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {
        sendString("countrange" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + minweight + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return numhits;
    }
    /** 
     * returns the sum of the hit weights on the specified alignment and chromosome between positions start and stop
     */
    public double getWeightRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException {
        return getWeightRange(alignid, chromid, start, stop, Float.NaN);
    }
    /** 
     * Returns the sum of the hit weights on the specified alignment and chromosome between positions start and stop.
     * Only includes hits whose weight is greater than minweight
     */
    public double getWeightRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {
        sendString("weightrange" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + minweight + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        return Double.parseDouble(readLine());
    }
    /** 
     * returns the sorted (ascending order) hit positions in the specified range of a chromosome,alignment pair.
     */ 
    public int[] getHitsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException {
        return getHitsRange(alignid, chromid, start,stop, Float.NaN);
    }
    /**
     * Returns the sorted hit positions whose weight is greater than minweight
     */
    public int[] getHitsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {
        sendString("gethitsrange" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + minweight + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        //        System.err.println("Going to read " + numhits + " hits");
        return Bits.readInts(numhits, instream, buffer, myorder);        
    }
    /**
     * Returns the hit weights in the specified range.  These will be in the same order
     * as the positions returned by getHitsRange();
     */
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException {
        return getWeightsRange(alignid, chromid, start,stop, Float.NaN);
    }
    /**
     * Returns weights in the range that are greater than minweight
     */
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {
        sendString("getweightsrange" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + minweight + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return Bits.readFloats(numhits, instream, buffer, myorder);        
    }
    /**
     * returns all hits on the specified chromosome
     */
    public int[] getHits(String alignid, String chromid)  throws IOException, ClientException {
        sendString("gethits" + "\n" + alignid + "\n" + chromid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return Bits.readInts(numhits, instream, buffer, myorder);        
    }
    /**
     * Returns all the weights on the specified chromosome
     */
    public float[] getWeights(String alignid, String chromid)  throws IOException, ClientException {
        sendString("getweights" + "\n" + alignid + "\n" + chromid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        int numhits = Integer.parseInt(readLine());
        return Bits.readFloats(numhits, instream, buffer, myorder);        
    }
    /**
     * returns a TreeMap from positions to counts representing a histogram
     * of the hits in a range with the specified binsize.  Bins with a count
     * of zero are not included in the output.
     * 
     * Ex: getHistgram("foo","1+",1,100,10)
     *     6: 5
     *    16: 0
     *    36: 30
     *
     * minweight is the minimum weight for reads to be included in the histogram.
     *
     * Normally, a read is only counted in a single bin as defined by its position (generally
     * the 5' end of the read).  A non-zero read-Extension counts the read in any
     * bin that you cross within that many bases of the read's position.  A negative value
     * goes backwards (smaller coordinates) and a larger value goes forward.  You need
     * to get the sign right depending on the strandedness of the chromosome that you're
     * working with.
     */
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {
        return getHistogram(alignid, chromid, start,stop,binsize, Float.NaN,0);
    }
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) 
        throws IOException, ClientException {
        sendString("histogram" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + binsize + "\n" + minweight + "\n" + readExtension + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            System.err.println("Asking for histogram said " + response);
            throw new ClientException(response);
        }
        int numints = Integer.parseInt(readLine());
        int out[] = Bits.readInts(numints, instream, buffer, myorder);
        TreeMap<Integer,Integer> output = new TreeMap<Integer,Integer>();
        for (int i = 0; i < out.length; i += 2) {
            output.put(out[i], out[i+1]);
        }
        return output;
    }
    /**
     * returns a histogram where the value is the summed read weights for the bin
     */
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {
        return getWeightHistogram(alignid, chromid, start,stop,binsize, Float.NaN,0);
    }
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) 
        throws IOException, ClientException {
        sendString("weighthistogram" + "\n" + alignid + "\n" + chromid + "\n" + start + "\n" + stop + "\n" + binsize + "\n" + minweight + "\n" + readExtension + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            System.err.println("Asking for weighthistogram said " + response);
            throw new ClientException(response);
        }
        int numints = Integer.parseInt(readLine());
        int pos[] = Bits.readInts(numints, instream, buffer, myorder);
        float weight[] = Bits.readFloats(numints, instream, buffer, myorder);
        TreeMap<Integer,Float> output = new TreeMap<Integer,Float>();
        for (int i = 0; i < pos.length; i ++) {
            output.put(pos[i], weight[i]);
        }
        return output;
    }
    /**
     * map from READ, WRITE, and ADMIN to lists of principals that have those privileges on the specified alignment
     */
    public Map<String,Set<String>> getACL(String alignid) throws IOException, ClientException {
        sendString("getacl\n" + alignid + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
        Map<String,Set<String>> output = new HashMap<String,Set<String>>();
        fillPartACL(output);
        fillPartACL(output);
        fillPartACL(output);
        return output;
    }
    /* fills one section of the acl output data structure.  A section
       is either read, write, or admin. 
    */
    private void fillPartACL(Map<String,Set<String>> output) throws IOException {
        String type = readLine();
        int entries = Integer.parseInt(readLine());
        Set<String> out = new HashSet<String>();
        while (entries-- > 0) {
            out.add(readLine());
        }
        output.put(type,out);
    }
    /**
     * Applies the specified ACLChangeEntry objects to the acl for this experiment/chromosome.
     */
    public void setACL(String alignid, Set<ACLChangeEntry> changes) throws IOException, ClientException {
        StringBuffer req = new StringBuffer("setacl\n" + alignid + "\n");
        for (ACLChangeEntry a : changes) {
            req.append(a.toString() + "\n");
        }    
        req.append("ENDREQUEST\n");
        sendString(req.toString());
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }   
    }
    public void addToGroup(String princ, String group) throws IOException, ClientException {
        sendString("addtogroup\n" + princ + "\n" + group + "\nENDREQUEST\n");
        String response = readLine();
        if (!response.equals("OK")) {
            throw new ClientException(response);
        }
    }
    /**
     * pack up and go home
     */
    public void close() {
        try {
            socket.setSoLinger(false,0);
            sendString("bye\nENDREQUEST\n");
            outstream.close();
            outstream = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            instream.close();
            instream = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
        try {
            socket.close();
            socket = null;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }    
    
}

/**
 * SASL callback handler for the authentication exchange
 */
class ClientCallbackHandler implements CallbackHandler {
    private String name, pass;
    public ClientCallbackHandler(String n, String p) {name = n; pass = p;}
    public void handle(Callback[] callbacks) {
        for (int i = 0; i < callbacks.length; i++) {
            if (callbacks[i] instanceof NameCallback) {
                NameCallback nc = (NameCallback)callbacks[i];
                nc.setName(name);
            }
            if (callbacks[i] instanceof PasswordCallback) {
                PasswordCallback pc = (PasswordCallback)callbacks[i];
                pc.setPassword(pass.toCharArray());
            }
        }
    }
}