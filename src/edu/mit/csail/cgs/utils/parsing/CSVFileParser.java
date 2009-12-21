/*
 * Created on Sep 25, 2007
 */
package edu.mit.csail.cgs.utils.parsing;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class CSVFileParser {

    private char quote, separator;
    private Pattern unquoter; 
    
    public CSVFileParser() { 
        separator = ',';
        quote = '"';
        unquoter = Pattern.compile(String.format("^%c(.*)%c$", quote, quote));
    }

    public CSVFileParser(char sep, char quo) { 
        separator = sep;
        quote = quo;
        unquoter = Pattern.compile(String.format("^%c(.*)%c$", quote, quote));
    }
    
    public String unquote(String f) { 
        Matcher m = unquoter.matcher(f);
        if(m.matches()) { 
            return m.group(1);
        } else { 
            return f;
        }
    }
    
    public Iterator<String[]> parseFile(File f) throws IOException { 
        return new CSVIterator(f);
    }
    
    public String[] parseLine(String line) {  
        LinkedList<Integer> separatorOffsets = new LinkedList<Integer>();
        
        boolean inQuote = false;
        for(int i = 0; i < line.length(); i++) { 
            char c = line.charAt(i);
            if(c == quote) { 
                inQuote = !inQuote;
            } else if (c==separator && !inQuote) { 
                separatorOffsets.addLast(i);
            }
        }
        
        return parseByOffsets(line, separatorOffsets);
    }
    
    private String[] parseByOffsets(String line, LinkedList<Integer> offs) { 
        int count = offs.size() + 1;
        String[] array = new String[count];
        
        int start = 0;
        int idx = 0;
        do { 
            int end = offs.isEmpty() ? line.length() : offs.removeFirst();
            String str = line.substring(start, end);
            array[idx++] = str;
            start = end + 1;
        } while(start < line.length());
        
        return array;
    }
    
    private class CSVIterator implements Iterator<String[]> { 
        
        private BufferedReader br;
        private String nextLine;
        
        public CSVIterator(File f) throws IOException { 
            br = new BufferedReader(new FileReader(f));
            nextLine = null;
            findNextLine();
        }
        
        public boolean hasNext() { return nextLine != null; }
        
        public String[] next() {
            String nl = nextLine;
            try {
                findNextLine();
            } catch (IOException e) {
                e.printStackTrace();
                nextLine = null;
                try {
                    br.close();
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
            }
            return parseLine(nl);
        }
        
        public void remove() { throw new UnsupportedOperationException(); }
        
        private void findNextLine() throws IOException { 
            do { 
                nextLine = br.readLine();
            } while(nextLine != null && nextLine.trim().length() == 0);
            if(nextLine != null) { 
                nextLine = nextLine.trim();
            } else { 
                br.close();
            }
        }
    }
}
