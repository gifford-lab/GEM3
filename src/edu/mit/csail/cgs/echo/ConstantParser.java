/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo;

import java.util.Map;
import java.util.regex.*;

import edu.mit.csail.cgs.ewok.types.SelfDescribingConstant;
import edu.mit.csail.cgs.ewok.types.ValueWrapper;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.species.*;

public interface ConstantParser<X extends SelfDescribingConstant> {
    public X parseConstant(String t);
    
    public static class RegionParser implements ConstantParser<ValueWrapper> {
    	
    	private Pattern p;
    	
    	public RegionParser() { 
    		p = Pattern.compile("(.+):(.+):(\\d+)-(\\d+)");
    	}

		public ValueWrapper parseConstant(String t) {
			Matcher m = p.matcher(t);
			if(m.matches()) { 
				String genomeName = m.group(1);
				try {
					Genome genome = Organism.findGenome(genomeName);
					String chrom = m.group(2);
					int start = Integer.parseInt(m.group(3));
					int end = Integer.parseInt(m.group(4));
					return new ValueWrapper(new Region(genome, chrom, start, end));
				} catch (NotFoundException e) {
					return null;
				}
			} else { 
				return null;
			}
		} 
    	
		public String toString() { 
			return "Region";
		}
    }
    
    public static class StringParser implements ConstantParser<ValueWrapper> {
        public ValueWrapper parseConstant(String t) {
            return new ValueWrapper(t);
        } 
        
        public String toString() { return "String"; }
    }
    
    public static class DoubleParser implements ConstantParser<ValueWrapper> {
        public ValueWrapper parseConstant(String t) {
            return new ValueWrapper(Double.parseDouble(t));
        } 
        
        public String toString() { return "Double"; }
    }
    
    public static class IntegerParser implements ConstantParser<ValueWrapper> {
        public ValueWrapper parseConstant(String t) {
            return new ValueWrapper(Integer.parseInt(t));
        } 
        
        public String toString() { return "Integer"; }
    }
    
    public static class GenomeParser implements ConstantParser<ValueWrapper> {
        public ValueWrapper parseConstant(String t) {
            try {
                Genome g = Organism.findGenome(t);
                return new ValueWrapper(g);
            } catch (NotFoundException e) {
                e.printStackTrace();
                return null;
            }
        }
        
        public String toString() { return "Genome"; }
    }
    
    public static class AnnotationParamsParser implements ConstantParser<ValueWrapper> {

        public String toString() { return "Annotation Params"; }
        
        public ValueWrapper parseConstant(String t) {
        	/*
            AnnotationParameters params = null;
            OptionsParser parser = new OptionsParser(t);
            Map<String,String> map = parser.getValueMap();
            int upstream = map.containsKey("up") ? Integer.parseInt(map.get("up")) : 8000;
            int downstream = map.containsKey("down") ? Integer.parseInt(map.get("down")) : 2000;
            int general = map.containsKey("general") ? Integer.parseInt(map.get("general")) : downstream;
            params = new RegionDistanceAnnotationParameters(upstream, downstream, general);
            return new ValueWrapper(params);
            */
        	throw new UnsupportedOperationException();
        } 
        
    }
}
