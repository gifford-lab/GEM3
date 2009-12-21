package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;

public class SimpleGeneMap implements GeneMap, Saveable {
    
    public static SimpleGeneMap loadSimpleGeneMap(File f) throws IOException { 
        DataInputStream dis = new DataInputStream(new FileInputStream(f));
        SimpleGeneMap sgm = new SimpleGeneMap(dis);
        dis.close();
        return sgm;
    }
	
	private Map<String,Set<String>> map;
	private Set<String> targets;
	private int size;
	
	public SimpleGeneMap() { 
        map = new HashMap<String,Set<String>>();
        targets = new HashSet<String>();
        size = 0;
	}
    
    public SimpleGeneMap(GeneMap base) { 
        map = new HashMap<String,Set<String>>();
        targets = new HashSet<String>();
        size = 0;
        
        for(String id : base.getDomainIDs()) { 
            Set<String> targets = base.mapID(id);
            for(String t : targets) { 
                addMapping(id, t);
            }
        }
    }
	
	public SimpleGeneMap(File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			if(line.length() > 0 && line.charAt(0) != '#') { 
				String[] array = line.split("\\s+");
				addMapping(array[0], array[1]);
			}
		}
		br.close();
	}
	
	public SimpleGeneMap(DataInputStream dis) throws IOException { 
		size = dis.readInt();
		map = new HashMap<String,Set<String>>();
		targets = new HashSet<String>();
		
		int s = dis.readInt();
		for(int i = 0; i < s; i++) { targets.add(dis.readUTF()); }
		
		s = dis.readInt();
		for(int i = 0; i < s; i++) { 
			String id = dis.readUTF();
			map.put(id, new HashSet<String>());
			int ss = dis.readInt();
			for(int j= 0; j < ss; j++) { 
				map.get(id).add(dis.readUTF());
			}
		}	
	}

	public void outputPairs(PrintStream ps) { 
		for(String k : map.keySet()) { 
			for(String v : map.get(k)) { 
				ps.println(k + "\t" + v);
			}
		}
	}
	

	public void save(DataOutputStream dos) throws IOException { 
		dos.writeInt(size);
		dos.writeInt(targets.size());
		for(String id : targets) { dos.writeUTF(id); }
		dos.writeInt(map.size());
		for(String id : map.keySet()) { 
			dos.writeUTF(id);
			dos.writeInt(map.get(id).size());
			for(String t : map.get(id)) { 
				dos.writeUTF(t);
			}
		}
	}
	
	public int size() { return size; }
	
	public void addMapping(String id, String target) { 
		if(!map.containsKey(id)) { map.put(id, new HashSet<String>()); }
		if(!map.get(id).contains(target)) { size += 1; }
		map.get(id).add(target);
		targets.add(target);
	}

	public Set<String> getDomainIDs() {
		return map.keySet(); 
	}

	public Set<String> getRangeIDs() {
		return targets;
	}

	public Set<String> mapID(String id) {
		return map.get(id);
	}

	public Set<String> mapIDs(Collection<String> ids) {
		HashSet<String> targs = new HashSet<String>();
		for(String id : ids) { targs.addAll(mapID(id)); }
		return targs;
	}

}
