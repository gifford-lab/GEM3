package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.RegionMapperFactory;

public class ByteSequenceGeneratorFactory implements RegionMapperFactory<Byte[]> {
    String type;

    public ByteSequenceGeneratorFactory() {
        type = "Byte[]";
    }
    public ByteSequenceGeneratorFactory(String t) {
        type = t;
    }

    public void setType(String t) {type = t;}
    public String getType() {return type;}
    public String getProduct() {return "Byte[]";}
    public Mapper<Region, Byte[]> getMapper(Genome g) {
        return getMapper(g,type);
    }

    public Mapper<Region, Byte[]> getMapper(Genome g, String type) {
        if (type == null) {
            throw new NullPointerException("Byte[] must have a type");
        } else {
            return new ByteSequenceGenerator(type);
        }
    }

    public Mapper<Region, Byte[]> getMapper(Genome g, String type, String k, String v) {
        if (type == null) {
            throw new NullPointerException("Byte[] must have a type");
        } else {
            return new ByteSequenceGenerator(type,k,v);
        }
    }
}
