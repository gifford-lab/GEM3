package edu.mit.csail.cgs.projects.readdb;

import org.apache.commons.cli.*;
import java.util.*;
import java.io.*;
import java.net.*;
import net.sf.samtools.*;
import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;



/**
 * Query hits from a remote BAM file.Reads chrom:start-stop:strand values from
 * stdin.  Repeats them on stdout along with the hit positions
 *
 * Usage:
 * Query --data http://nanog.csail.mit.edu/readdb-test/foo.bam --index foo.index [--quiet]
 *
 * --quiet means don't print any output.  This is useful for testing the query performance
 * without worrying about the time it takes to print the output.
 *
 * --data and --index can be either local or remote
 *
 */

public class BAMQuery {

    private String data, index;
    private SAMFileReader reader;
    private boolean quiet, weights, noheader;
    private int histogram;
    public static void main(String args[]) throws Exception {
        BAMQuery query = new BAMQuery();
        query.parseArgs(args);
        query.run(System.in);

    }
    public void parseArgs(String args[]) throws IllegalArgumentException, ParseException {
        Options options = new Options();
        options.addOption("q","quiet",false,"quiet: don't print output");
        options.addOption("d","data",true,"url for data");
        options.addOption("i","index",true,"url for index file");
        options.addOption("w","weights",false,"get and print weights in addition to positions");
        options.addOption("H","histogram",true,"produce a histogram with this binsize instead of printing all read positions");
        options.addOption("N","noheader",false,"skip printing the query header");
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse( options, args, false );            
        quiet = line.hasOption("quiet");
        weights = line.hasOption("weights");
        noheader = line.hasOption("noheader");
        if (line.hasOption("histogram")) {
            histogram = Integer.parseInt(line.getOptionValue("histogram"));
        } else {
            histogram = 0;
        }
        if (line.hasOption("data")) {
            data = line.getOptionValue("data");
        } else {
            throw new IllegalArgumentException("Must provide --data");
        }
        if (line.hasOption("index")) {
            index = line.getOptionValue("index");
        } else {
            throw new IllegalArgumentException("Must provide --index");
        }
    }
    public BAMQuery() {}
    public SAMFileReader createReader() throws IOException, URISyntaxException {
        URI indexURI = new URI(index);
        String indexFilename = null;
        if (indexURI.getScheme().equals("file")) {
            indexFilename = indexURI.getPath();
        } else if (indexURI.getScheme().equals("http")) {
            DefaultHttpClient httpclient = new DefaultHttpClient();
            HttpGet httpget = new HttpGet(indexURI);
            HttpResponse response = httpclient.execute(httpget);
            HttpEntity entity = response.getEntity();

            if (entity != null) {
                    entity.consumeContent();
                    FileOutputStream os = new FileOutputStream(indexFilename);
                    entity.writeTo(os);
                    os.close();
            }
            httpclient.getConnectionManager().shutdown();        
        } else {
            indexFilename = index; // hope for the best here
        }        
        return new SAMFileReader(new URL(data), new File(indexFilename), false);       
    }
    public void run(InputStream instream) throws IOException, URISyntaxException {
        SAMFileReader bam = createReader();
        //        bam.enableIndexCaching(true);
        BufferedReader reader = new BufferedReader(new InputStreamReader(instream));
        String line = null;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("[\\:]");
            String chr = pieces[0];
            Boolean strand = pieces.length >= 3 ? pieces[2].equals("-") : null;
            pieces = pieces[1].split("\\-");
            int start = Integer.parseInt(pieces[0]);
            int stop = Integer.parseInt(pieces[1]);
            
            SAMRecordIterator iter = bam.query(chr, start, stop, false);

            if (histogram > 0) {
                int hist[] = new int[(stop - start) / histogram + 1];
                while (iter.hasNext()) {
                    SAMRecord record = iter.next();
                    if (strand != null && strand != record.getReadNegativeStrandFlag()) {
                        continue;
                    }
                    int pos = record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart();
                    int bin = (pos - start) / histogram;
                    if (bin < hist.length) {
                        hist[bin]++;
                    }
                }
                if (!quiet) {
                    for (int i = 0; i < hist.length; i++) {
                        System.out.println(String.format("%d\t%d", start + i*histogram + histogram/2, hist[i]));
                    }
                }
            } else {
                while (iter.hasNext()) {
                    SAMRecord record = iter.next();
                    if (strand != null && strand != record.getReadNegativeStrandFlag()) {
                        continue;
                    }
                    if (!quiet) {
                        System.out.println(String.format("chrom %d, pos %d, %s, len %d",
                                                         record.getReferenceName(),
                                                         record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart(),
                                                         record.getReadNegativeStrandFlag() ? "-" : "+",
                                                         record.getReadLength()));
                    }


                }

            }
        }
        reader.close();
        bam.close();
    }






}