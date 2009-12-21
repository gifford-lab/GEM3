package edu.mit.csail.cgs.ewok.verbs.sequence;

import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class StringToLongMapper implements Mapper<String,Long> {

    public Long execute(String a) {
        return SequenceUtils.StringToLong(a);
    } 	   
}