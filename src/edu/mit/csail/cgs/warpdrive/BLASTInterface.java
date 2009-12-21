/*
 * Created on Mar 1, 2007
 */
package edu.mit.csail.cgs.warpdrive;

import java.util.*;
import java.io.*;
import java.net.*;

import javax.xml.parsers.SAXParser;

import org.xml.sax.*;
import org.xml.sax.helpers.*;

import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.xml.XMLTreeElement;

public class BLASTInterface {
    
    public static void main(String[] args) { 
        try {
            String genomeName = args[0];
            String urlString = args[1];
            BLASTInterface blaster = new BLASTInterface(genomeName, urlString);
            String queryResponse = blaster.doQuery(args[2]);
            
            BLASTHitHandler handler = new BLASTHitHandler();
            
            XMLReader oParser = XMLReaderFactory.createXMLReader(); 
            
            oParser.setContentHandler(handler);
            
            StringReader reader = new StringReader(queryResponse);
            oParser.parse(new InputSource(reader));
            
            System.out.println("Hits:");
            Genome genome = Organism.findGenome(genomeName);
            for(ScoredStrandedRegion hit : handler.extractHitRegions(genome)) { 
                System.out.println(hit);
            }

        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    private String genomeName;
    private URL blastURL;
    private String agent, type;
    private String encType = "multipart/form-data";
    
    public BLASTInterface(String genome, String urlString) throws MalformedURLException {
        genomeName = genome;
        blastURL = new URL(urlString);
        agent = "Java agent";
    }
    
    public String getGenomeName() { return genomeName; }
    
    private String buildEncodedData(String query) { 
        StringBuilder sb = new StringBuilder();
        sb.append("PROGRAM=blastn");
        sb.append("&");
        sb.append("DATALIB=" + genomeName);
        sb.append("&");
        sb.append("SEQUENCE="); sb.append(query);
        sb.append("&");
        sb.append("ALIGNMENT_VIEW=7");
        return sb.toString();
    }
    
    public String doQuery(String query) throws IOException {
        String encData = buildEncodedData(query);
        
        HttpURLConnection cxn = (HttpURLConnection)blastURL.openConnection();
        cxn.setRequestMethod("POST");
        cxn.setRequestProperty("User-Agent", agent);
        cxn.setRequestProperty("Content-Type", encType);
        cxn.setRequestProperty("Content-Length", String.valueOf(encData.length()));
        cxn.setDoOutput(true);

        OutputStream os = cxn.getOutputStream();
        os.write( encData.getBytes() );
    
        int rc = cxn.getResponseCode();
        System.out.println("Response: " + rc);
        
        InputStream is = cxn.getInputStream();
        BufferedReader br = new BufferedReader(new InputStreamReader(is));
        
        StringBuilder doc = new StringBuilder();
        
        String line = null;
        int i = 0;
        while((line = br.readLine()) != null) { 
            doc.append(line);
            i++;
        }
        
        br.close();
        cxn.disconnect();
        
        return doc.toString();
    }
}
