/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jsplice.exception.Log;
import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class PatternMotif implements Comparable<PatternMotif>{

  public String pattern;
  public int quantityRef;
  public int quantityAlt;
  public int quantityBen;
  public int quantityPat;
  public int shift;

  /**
   * 
   */
  public PatternMotif(String patternP, int quantityAbsP, int quantityConP, int shiftP) {
    this.pattern = patternP;
    this.quantityRef = quantityAbsP;
    this.quantityAlt = quantityConP;
    this.shift = shiftP;
  }

  public PatternMotif(PatternMotif patternP) {
    this.pattern = patternP.pattern;
    this.quantityRef = patternP.quantityRef;
    this.quantityAlt = patternP.quantityAlt;
    this.quantityBen = patternP.quantityBen;
    this.quantityPat = patternP.quantityPat;
    this.shift = patternP.shift;
  }
  /* (non-Javadoc)
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return pattern + " " + quantityBen + " " + getQuantityRelative();
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals()
   */
  @Override
  public boolean equals(Object patternP) {
    if (patternP == null) {
      return false;
    } else if (patternP == this) {
      return true;
    }
    else if ((patternP instanceof PatternMotif)) {
      PatternMotif patternMotif = (PatternMotif) patternP;
      return this.pattern.equals(patternMotif.pattern);
    }
//    else if ((patternP instanceof String)) {
//      String str = (String) patternP;
//      return this.pattern.equals(str);
//    }
    return super.equals(patternP);
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
    return (double) quantityRef / (quantityAlt + 1);
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

  public static int align(String str1, String str2){
    int matchesMax = 0;
    int posMax = -1;
//    String str1 = pattern1.pattern;
//    String str2 = pattern2.pattern;
    for (int pos = str1.length()-1; pos > -str2.length(); pos--) {
      int from1 = pos > 0 ? pos : 0;
      int to1 = pos + str2.length() < str1.length() ? pos + str2.length() : str1.length();
      int from2 = -pos > 0 ? -pos : 0;
      int to2 = -pos + str1.length() < str2.length() ? -pos + str1.length() : str2.length();
//      System.out.println(from1 + " > " + to1 + " of " + str1);
//      System.out.println(from2 + " > " + to2 + " of " + str2);
      int matches = Functions.countMatchingChar(str1.substring(from1, to1), str2.substring(from2, to2));
      //			System.out.println(matches);
      if (matches > matchesMax) {
        matchesMax = matches;
        posMax = pos;
      } else if (matches == matchesMax) {
        Log.add(posMax + " and " + pos + " are possible possitions.", 2);
      }
    }
    return posMax;
  }
}
