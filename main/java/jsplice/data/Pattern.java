/**
 * 
 */
package jsplice.data;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Pattern implements Comparable<Pattern>{

	public String pattern;
	public int quantityAbs;
	public int quantityCon;
	public int quantityUnique;
	
	/**
	 * 
	 */
	public Pattern(String patternP, int quantityAbsP, int quantityConP) {
		this.pattern = patternP;
		this.quantityAbs = quantityAbsP;
		this.quantityCon = quantityConP;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return pattern + " " + quantityUnique + " " + getQuantityRelative();
	}
	
	public boolean equals(Pattern patternP) {
	    if (patternP == null) 
	    	return false;
	    if (patternP == this) 
	    	return true;
	    if (!(patternP instanceof Pattern))
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
	public int compareTo(Pattern patternP) {
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
	public boolean contains(Pattern patternP) {
		return pattern.contains(patternP.pattern);
	}

	/**
	 * @param patternP
	 * @return
	 */
	public boolean contains(String patternP) {
		return pattern.contains(patternP);
	}
}
