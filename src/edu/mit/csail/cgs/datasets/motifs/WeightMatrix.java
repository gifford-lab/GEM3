package edu.mit.csail.cgs.datasets.motifs;

import java.util.*;
import java.sql.*;
import java.text.DecimalFormat;
import java.text.FieldPosition;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.NotFoundException;

/* weight matrix representation as float[][]:
   the first index is the position in the matrix.
   the second index is the letter.  The size of the second dimension
   is MAXLETTERVAL, but only eight of the indexes are valid: A, C, T, G, and the
   four corresponding lower case letters.  The lower case entries should have the same value
   as the upper case entries.
   This uses more space but makes for quick lookups because you can
   use char values as the index */

public class WeightMatrix {

    public static int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
                                              Math.max(Math.max('a','c'),Math.max('t','g')))
        + 1;
    public static char[] allLetters = {'A','a','C','c','T','t','G','g'};
    public static char[] letters = {'A','C','T','G'};
    public static char[] revCompLetters = {'T','G','A','C'};

    public float[][] matrix;
    public String name, version, type, species;
    public int dbid, speciesid;
    public boolean hasdbid, hasspeciesid; // set to true iff dbid is valid
    public boolean islogodds;
    public int bgModelID = -1;
    
    WeightMatrix(ResultSet wmData, ResultSet wmColData) throws SQLException {  
    	dbid = wmData.getInt(1);
    	speciesid = wmData.getInt(2);
    	name = wmData.getString(3);
    	version = wmData.getString(4);
    	type = wmData.getString(5);
    	bgModelID = wmData.getInt(6);
    	if (wmData.wasNull()) {
    	  bgModelID = -1;
    	}
    	
    	hasdbid = true; hasspeciesid = true;
    	islogodds = false;
    	
    	int maxPos = -1;
    	Map<Integer,Map<Character,Float>> weights = new HashMap<Integer,Map<Character,Float>>();
    	
    	while(wmColData.next()) { 
    		int pos = wmColData.getInt(1);
    		char posChar = wmColData.getString(2).charAt(0);
    		float weight = wmColData.getFloat(3);
            if(weight < 0.0) { islogodds = true; }

            if(!weights.containsKey(pos)) { 
    			maxPos = Math.max(maxPos, pos);
    			weights.put(pos, new HashMap<Character,Float>());
    		}
    		
    		weights.get(pos).put(posChar, weight);
    	}

        float nullWeight = islogodds ? (float)-9999.0 : (float)0.0;
    	
    	matrix = new float[maxPos+1][MAXLETTERVAL];
        for(int i = 0; i < matrix.length; i++) { 
            for(int j = 0; j < matrix[j].length; j++) { 
                matrix[i][j] = nullWeight;
            }
        }
		for(int j = 0; j < letters.length; j++) { 
			char letter = letters[j];
			for(int i = 0; i < matrix.length; i++) { 
				if(weights.containsKey(i) && weights.get(i).containsKey(letter)) { 
					matrix[i][Character.toLowerCase(letter)] = weights.get(i).get(letter);
					matrix[i][Character.toUpperCase(letter)] = weights.get(i).get(letter);
				} else { 
					matrix[i][Character.toLowerCase(letter)] = nullWeight;
					matrix[i][Character.toUpperCase(letter)] = nullWeight;
				}
    		}
    	}
    }    

    public WeightMatrix (int length) {
        dbid = -1;
        hasdbid = false;
        hasspeciesid = false;
        matrix = new float[length][MAXLETTERVAL];                
    }

    public static Collection<WeightMatrix> getAllWeightMatrices() {
        try {
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select m.id, m.species, m.name, m.version, m.type, m.bg_model_map_id, c.position, c.letter, c.weight from "
                                                        + "weightmatrix m, weightmatrixcols c where "
                                                        + " m.id = c.weightmatrix order by c.weightmatrix, c.position desc");
            ResultSet rs = ps.executeQuery();

            Collection<WeightMatrix> matrices = WeightMatrix.getWeightMatrices(rs);

            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return matrices;
        } catch (SQLException ex){ 
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    /**
     * Gets all the matrices specified in the result set which is
     * the result of joining weightmatrix and weightmatrixcols, eg
     * select m.id, m.species, m.name, m.version, m.type, m.bg_model_map_id, c.position, c.letter, c.weight from weightmatrix m, 
     * weightmatrixcols c where m.id = c.weightmatrix order by c.weightmatrix, c.position desc
     *
     * sorting by descending position is critical so that the first row of a new matrix gives
     * the largest index (ie, the length of the matrix)
     */
    public static Collection<WeightMatrix> getWeightMatrices(ResultSet rs) throws SQLException {
        Map<Integer,WeightMatrix> output = new HashMap<Integer,WeightMatrix>();
        while (rs.next()) {
            int id = rs.getInt(1);
            if (!output.containsKey(id)) {
                WeightMatrix m = new WeightMatrix(rs.getInt(7) + 1);
                m.dbid = rs.getInt(1);
                m.hasdbid = true;
                m.speciesid = rs.getInt(2);
                if (m.speciesid > 0) {
                  m.hasspeciesid = true;
                }
                m.name = rs.getString(3);
                m.version = rs.getString(4);
                m.type = rs.getString(5);
                m.bgModelID = rs.getInt(6);
                if ((m.bgModelID == 0) && rs.wasNull()) {
                  m.bgModelID = -1;
                }
                output.put(id,m);
            }
            output.get(id).matrix[rs.getInt(7)][rs.getString(8).charAt(0)] = rs.getFloat(9);
        }
        System.err.println("Returning " + output.size() + " from WeightMatrix.getWeightMatrices(ResultSet)");
        return output.values();
    }

    /* creates a new weightMatrixObject based on the provided database identifier */
    public static WeightMatrix getWeightMatrix(int dbid) throws NotFoundException {
        WeightMatrix matrix = null;
        try {
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select max(position) from weightmatrixcols where weightmatrix = ?");
            ps.setInt(1,dbid);
            int length;
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                length = rs.getInt(1) + 1;
            } else {
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
                throw new NotFoundException("Can't find WMID " + dbid);
            }
            rs.close();
            ps.close();
            matrix = new WeightMatrix(length);
            matrix.dbid = dbid;
            matrix.hasdbid = true;
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < MAXLETTERVAL; j++) {
                    matrix.matrix[i][j] = Float.NaN;
                }
            }
            ps = cxn.prepareStatement("select position, letter, weight from weightmatrixcols where weightmatrix = ? order by position");
            ps.setInt(1,dbid);
            rs=ps.executeQuery();
            while (rs.next()) {
              float weight = rs.getFloat(3);
              matrix.matrix[rs.getInt(1)][rs.getString(2).charAt(0)] = weight;
              matrix.matrix[rs.getInt(1)][Character.toLowerCase(rs.getString(2).charAt(0))] = weight;
                
              if(weight < 0.0) { matrix.islogodds = true; }
            }
            rs.close();
            ps.close();
            ps = cxn.prepareStatement("select name, version, type, species, bg_model_map_id from weightmatrix where id = ?");
            ps.setInt(1,dbid);
            rs = ps.executeQuery();
            rs.next();
            matrix.name = rs.getString(1);
            matrix.version = rs.getString(2);
            matrix.type = rs.getString(3);
            matrix.speciesid = rs.getInt(4);
            matrix.hasspeciesid = true;
            matrix.bgModelID = rs.getInt(5);
            if ((matrix.bgModelID == 0) && rs.wasNull()) {
              matrix.bgModelID = -1;
            }
            rs.close();
            ps.close();            
            DatabaseFactory.freeConnection(cxn);     
        } catch (SQLException ex) {
            throw new NotFoundException("Can't find WM " + dbid,ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Can't connect to annotations datasource",ex);
        }
        return matrix;
    }

    public int length () {return matrix.length;}
    public String getName(){return name;}
    public String getVersion(){return version;}
    /* returns the maximum possible score that a sequence
       could have against the specified matrix */
    public double getMaxScore() {
        float maxscore = 0;
        for (int i = 0; i < matrix.length; i++) {
            float max = Float.NEGATIVE_INFINITY;
            if (matrix[i]['A'] > max) {
                max = matrix[i]['A'];
            }
            if (matrix[i]['C'] > max) {
                max = matrix[i]['C'];
            }
            if (matrix[i]['G'] > max) {
                max = matrix[i]['G'];
            }
            if (matrix[i]['T'] > max) {
                max = matrix[i]['T'];
            }

            maxscore += max;
        }
        return maxscore;
    }
    public double getMinScore() {
        float minscore = 0;
        for (int i = 0; i < matrix.length; i++) {
            float min = Float.POSITIVE_INFINITY;
            if (matrix[i]['A'] < min) {
                min = matrix[i]['A'];
            }
            if (matrix[i]['C'] < min) {
                min = matrix[i]['C'];
            }
            if (matrix[i]['G'] < min) {
                min = matrix[i]['G'];
            }
            if (matrix[i]['T'] < min) {
                min = matrix[i]['T'];
            }

            minscore += min;
        }
        return minscore;
    }
    public boolean equals(Object other) {
        if (other instanceof WeightMatrix) {
            WeightMatrix o = (WeightMatrix)other;
            if (o.hasdbid && this.hasdbid) {
                return o.dbid == this.dbid;
            } else {
                return (o.name != null &&
                        this.name != null &&
                        this.name.equals(o.name)
                        && o.version != null
                        && this.version != null
                        && this.version.equals(o.version));
            }
        } else {
            return false;
        }
    }
    public int hashCode() {
        if (hasdbid) {
            return dbid;
        }
        return (name.hashCode() + version.hashCode());
    }
    public int getDBID() {return dbid;}
    public String toString() {
        return name + ", " + version + ", " + type;
    }

    public WeightMatrix subMatrix(int start, int length) {
        WeightMatrix out = new WeightMatrix(length);
        for (int i = start; i < start + length; i++) {
            out.matrix[i - start] = matrix[i];
        }
        out.name = name + "(" + start + "-" + (start +length + 0.0) + ")";
        out.version = version;
        out.type = type;
        return out;
    }

    /* examines the weight matrix to determine if it's a log-odds matrix or not */
    public boolean setLogOdds() {
        islogodds = false;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < allLetters.length; j++) {
                if (matrix[i][allLetters[j]] < 0) {
                    islogodds = true;
                    return true;
                }
            }
        }
        return false;
    }

    /* converts this matrix to log-odds form */
    public void toLogOdds() {
        setLogOdds();
        if (islogodds) {return;}
        islogodds = true;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < allLetters.length; j++) {
              MarkovBackgroundModel bgModel = null;
            	
            	if (bgModelID != -1) {
            	  try {
            	    bgModel = BackgroundModelLoader.getMarkovModel(bgModelID);
            	  }
            	  catch (NotFoundException nfex) {
            	    nfex.printStackTrace();
            	  }
            	  catch (SQLException sqlex) {
            	    sqlex.printStackTrace();
            	  }                
            	}
            	
            	if (bgModel != null) {
            	  matrix[i][allLetters[j]] = (float)Math.log(Math.max(matrix[i][allLetters[j]], .000001) / bgModel.getMarkovProb(("" + allLetters[j]).toUpperCase()));
            	}
            	else {
            		matrix[i][allLetters[j]] = (float)Math.log(Math.max(matrix[i][allLetters[j]], .000001) / .25);
            	}
            }
        }        
    }
    /* converts this matrix to frequencies.
       This conversion is in-exact because it uses a default background
       model */
    public void toFrequency() {
        setLogOdds();
        if (!islogodds) {return;}
        islogodds = false;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < allLetters.length; j++) {
                matrix[i][allLetters[j]] = (float)(Math.exp(matrix[i][allLetters[j]]) * .25);
            }
        }
    }
    /**
     * Modifies this WeightMatrix by normalizing the entries.  The matrix
     * must be a frequency matrix.  The result is that the frequencies at each
     * position in the matrix sum to 1.  This method is useful if you've read
     * in a matrix of counts and need to convert to frequencies.
     */
    public void normalizeFrequencies() {
        setLogOdds();
        if (islogodds) {
            return;
        }
        for (int position = 0; position < matrix.length; position++) {
            double sum = 0;
            for (int j = 0; j < letters.length; j++) {
                sum += matrix[position][letters[j]];
            }
            for (int j = 0; j < letters.length; j++) {
                matrix[position][letters[j]] = (float)(matrix[position][letters[j]] / sum);
            }
        }
    }

    /* looks up the database id for the specified weight matrix */
    public static int getWeightMatrixID(int speciesid,
                                        String wmname,
                                        String wmversion) throws NotFoundException {
        int wmid = -1;
        try {
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select id from weightmatrix where species = ? " + 
                                                        "and name = ? and version = ?");
            ps.setInt(1,speciesid);
            ps.setString(2,wmname);
            ps.setString(3,wmversion);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                wmid = rs.getInt(1);
            } else {
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);     
                throw new NotFoundException("Can't find WM " + wmname + ", " + wmversion + " in species " + 
                                            speciesid);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);     
        } catch (SQLException ex) {
            throw new NotFoundException("Can't find WM " + wmname + ", " + wmversion + " in species " + 
                                        speciesid,ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Can't connect to annotations datasource",ex);
        }
        return wmid;
    }

    /* returns a string showing the full matrx in all its numeric glory */
    public static String printMatrix(WeightMatrix matrix) {
        int length = matrix.length();
        StringBuffer out = new StringBuffer();
        DecimalFormat format = new DecimalFormat("##.##");
        out.append("\t");
        for (int i = 0; i < length;i++) {
            out.append(i + "\t");
        }
        out.append("\n");
        out.append("A\t");
        for (int i = 0; i < length; i++) {            
            out.append(format.format(matrix.matrix[i]['A']) + "\t");
        }
        out.append("\n");
        out.append("C\t");
        for (int i = 0; i < length; i++) {
            out.append(format.format(matrix.matrix[i]['C']) + "\t");
        }
        out.append("\n");
        out.append("G\t");
        for (int i = 0; i < length; i++) {
            out.append(format.format(matrix.matrix[i]['G']) + "\t");
        }
        out.append("\n");
        out.append("T\t");
        for (int i = 0; i < length; i++) {
            out.append(format.format(matrix.matrix[i]['T']) + "\t");
        }
        out.append("\n");
        return out.toString();
    }

    /* returns a string showing a text representation of the motif */
    public static String printMatrixLetters(WeightMatrix matrix) {
        char[][] out = new char[4][matrix.length()];
        for (int i = 0; i < out.length; i++) {
            for (int j = 0; j < out[0].length; j++) {
                out[i][j] = ' ';
            }
        }
        Character letters[] = {'A','C','G','T'};
        WMLetterCmp cmp = new WMLetterCmp(matrix);
        double minval = matrix.setLogOdds() ? 0 : .25;
        for (int i = 0; i < matrix.length(); i++) {
            cmp.setIndex(i);
            Arrays.sort(letters,cmp);
            for (int j = 0; j < 4; j++) {
                if (matrix.matrix[i][letters[j]] > minval) {
                    out[j][i] = letters[j];
                } else {
                    break;
                }
            }
        }
        return new String(out[0]) + "\n" +
            new String(out[1]) + "\n" +
            new String(out[2]) + "\n" +
            new String(out[3]);
    }
    /* returns the four letters ordered by their LL at the specified index */
    public static Character[] getLetterOrder(WeightMatrix wm, int index) {
        Character letters[] = {'A','C','G','T'};
        WMLetterCmp cmp = new WMLetterCmp(wm);
        cmp.setIndex(index);
        Arrays.sort(letters,cmp);
        return letters;
    }

}
class WMLetterCmp implements Comparator<Character> {
    WeightMatrix matrix;
    int index;
    public WMLetterCmp(WeightMatrix m) {matrix = m;}
    public void setIndex(int i) {index = i;}
    public int compare(Character a, Character b) {
        return Double.compare(matrix.matrix[index][b],matrix.matrix[index][a]);
    }
}
