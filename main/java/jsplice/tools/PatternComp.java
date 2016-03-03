/**
 * 
 */
package jsplice.tools;

import java.util.Comparator;

import jsplice.data.PatternMotif;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class PatternComp {
  public static Comparator<PatternMotif> desc(final Comparator<PatternMotif> sort) {
    return new Comparator<PatternMotif>() {
      public int compare(PatternMotif c1, PatternMotif c2) {
        return -1 * sort.compare(c1, c2);
      }
    };
  }
  
  public static Comparator<PatternMotif> order(final Comparator<PatternMotif> sort1, final Comparator<PatternMotif> sort2) {
    return new Comparator<PatternMotif>() {
      public int compare(PatternMotif o1, PatternMotif o2) {
        int res1 = sort1.compare(o1, o2);
        int res2 = sort2.compare(o1, o2);
        return res1 != 0 ? res1 : res2;
      }
    };
  }
}

enum AttrPat implements Comparator<PatternMotif> {
  QTY_REF {
    public int compare(PatternMotif p1, PatternMotif p2) {
      return Double.compare(p1.getQuantityRef(), p2.getQuantityRef());
    }
  },
  QTY_BEN {
    public int compare(PatternMotif p1, PatternMotif p2) {
      return Double.compare(p1.getQuantityBen(), p2.getQuantityBen());
    }
  },
  QTY_REF_REL {
    public int compare(PatternMotif p1, PatternMotif p2) {
      return Double.compare(p1.getQuantityRefRelative(), p2.getQuantityRefRelative());
    }
  },
  QTY_BEN_REL {
    public int compare(PatternMotif p1, PatternMotif p2) {
      return Double.compare(p1.getQuantityBenRelative(), p2.getQuantityBenRelative());
    }
  },
  LENGTH {
    public int compare(PatternMotif p1, PatternMotif p2) {
      return Integer.compare(p1.length(), p2.length());
    }
  };
}
