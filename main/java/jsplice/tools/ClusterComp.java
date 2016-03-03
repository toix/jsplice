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
public class ClusterComp {
  public static Comparator<Cluster> desc(final Comparator<Cluster> sort) {
    return new Comparator<Cluster>() {
      public int compare(Cluster c1, Cluster c2) {
        return -1 * sort.compare(c1, c2);
      }
    };
  }
  public static Comparator<Cluster> order(final Comparator<Cluster> sort1, final Comparator<Cluster> sort2) {
    return new Comparator<Cluster>() {
      public int compare(Cluster o1, Cluster o2) {
        int res1 = sort1.compare(o1, o2);
        int res2 = sort2.compare(o1, o2);
        return res1 != 0 ? res1 : res2;
      }
    };
  }
}


enum AttrCl implements Comparator<Cluster> {
  QTY_REF {
    public int compare(Cluster c1, Cluster c2) {
      return Double.compare(c1.getPatternCore().getQuantityRef(), c2.getPatternCore().getQuantityRef());
    }
  },
  QTY_REF_REL {
    public int compare(Cluster c1, Cluster c2) {
      return Double.compare(c1.getPatternCore().getQuantityRefRelative(), c2.getPatternCore().getQuantityRefRelative());
    }
  },
  QTY_BEN {
    public int compare(Cluster c1, Cluster c2) {
      return Double.compare(c1.getQuantityBen(), c2.getQuantityBen());
    }
  },
  QTY_BEN_REL {
    public int compare(Cluster c1, Cluster c2) {
      return Double.compare(c1.getQuantityBenRelative(), c2.getQuantityBenRelative());
    }
  },
  LENGTH {
    public int compare(Cluster c1, Cluster c2) {
      return Double.compare(c1.getPatternCore().length(), c2.getPatternCore().length());
    }
  };
}
