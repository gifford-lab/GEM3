/*
 * Created on Sep 7", 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.BIND;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class BindReference {
    
    private int pmid;
    private String citation;
    private int method;
    
    public BindReference(int pmid, int method) {
        this.pmid = pmid;
        citation = null;
        this.method =method;
        if(method < -1 || method >= Methods.length) { 
            throw new IllegalArgumentException(String.valueOf(method));
        }
    }
    
    public BindReference(int pmid, String citation) { 
        this.pmid = pmid;
        this.citation = citation;
        method = -1;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof BindReference)) { return false; }
        BindReference br = (BindReference)o;
        if(pmid != br.pmid) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17 + pmid; 
        code *= 37;
        return code;
    }
    
    public String toString() { 
        String str = "[" + pmid;
        if(citation != null) { str += " (" + citation + ")"; }
        if(method != -1) { str += " <" + getMethod() + ">"; }
        str += "]";
        return str;
    }
    
    public int getPMID() { return pmid; }
    public String getCitation() { return citation; }
    public int getMethodID() { return method; }
    public String getMethod() { return Methods[method]; }

    public static final String[] Methods = {
        "not-specified",
        "alanine-scanning",
        "affinity-chromatography",
        "atomic-force-microscopy",
        "autoradiography",
        "competition-binding",
        "cross-linking",
        "deuterium-hydrogen-exchange",
        "electron-microscopy",
        "electron-spin-resonance",
        "elisa",
        "equilibrium-dialysis",
        "fluorescence-anisotropy",
        "footprinting",
        "gel-retardation-assays",
        "gel-filtration-chromatography",
        "hybridization",
        "immunoblotting",
        "immunoprecipitation",
        "immunostaining",
        "interaction-adhesion-assay",
        "light-scattering",
        "mass-spectrometry",
        "membrane-filtration",
        "monoclonal-antibody-blockade",
        "nuclear-translocation-assay",
        "phage-display",
        "reconstitution",
        "resonance-energy-transfer",
        "site-directed-mutagenesis",
        "sucrose-gradient-sedimentation",
        "surface-plasmon-resonance-chip",
        "transient-coexpression",
        "three-dimensional-structure",
        "two-hybrid-test",
        "allele-specific-complementation",
        "far-western",
        "colocalization",
        "other" };
    public static Map<String,Integer> methodName2Index;
    static { 
        methodName2Index = new HashMap<String,Integer>();
        for(int i = 0; i < Methods.length; i++) { 
            methodName2Index.put(Methods[i], i);
        }
    }
}
