/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class PatternMotif implements Comparable<PatternMotif>{

	public String pattern;
	public int quantityAbs;
	public int quantityCon;
	public int quantityUnique;
	
	/**
	 * 
	 */
	public PatternMotif(String patternP, int quantityAbsP, int quantityConP) {
		this.pattern = patternP;
		this.quantityAbs = quantityAbsP;
		this.quantityCon = quantityConP;
	}
	
	public PatternMotif(PatternMotif patternP) {
		this.pattern = patternP.pattern;
		this.quantityAbs = patternP.quantityAbs;
		this.quantityCon = patternP.quantityCon;
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return pattern + " " + quantityUnique + " " + getQuantityRelative();
	}
	
	public boolean equals(PatternMotif patternP) {
	    if (patternP == null) 
	    	return false;
	    if (patternP == this) 
	    	return true;
	    if (!(patternP instanceof PatternMotif))
	    	return false;
    	return this.pattern.equals(patternP.pattern);
	}
	
	public int length(){
		return pattern.length();
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public int compareTo(PatternMotif patternP) {
		int lengthDelta = Integer.compare(patternP.length(), length());
		int quantytyRelDelta = Double.compare(patternP.getQuantityRelative(), getQuantityRelative());
		return lengthDelta != 0 ? lengthDelta : quantytyRelDelta;
	}
	
	public double getQuantityRelative() {
		return (quantityAbs - 1) / (quantityCon + 1);
	}

	/**
	 * @param patternP
	 * @return
	 */
	public boolean contains(PatternMotif patternP) {
		return pattern.contains(patternP.pattern);
	}

	/**
	 * @param patternP
	 * @return
	 */
	public boolean contains(String patternP) {
		return pattern.contains(patternP);
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#clone()
	 */
	@Override
	public PatternMotif clone() {
		return new PatternMotif(this);
	}

	public static ArrayList<PatternMotif> clone(List<PatternMotif> list) {
		ArrayList<PatternMotif> clone = new ArrayList<PatternMotif>(list.size());
	    for(PatternMotif item: list) {
	    	clone.add(new PatternMotif(item));
	    }
	    return clone;
	}

	/**
	 * @param list
	 * @return
	 */
	public static HashMap<PatternMotif, PatternMotif> cloneToMap(ArrayList<PatternMotif> list) {
		HashMap<PatternMotif, PatternMotif> clone = new HashMap<PatternMotif, PatternMotif>(list.size());
	    for(PatternMotif item: list) {
	    	clone.put(item, new PatternMotif(item));
	    }
	    return clone;
	}
}
