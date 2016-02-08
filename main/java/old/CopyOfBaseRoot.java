/**
 * 
 */
package old;

import java.util.ArrayList;

import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class CopyOfBaseRoot {

	public String bases;
	public int quantity = 0;
	public ArrayList<Base> subBasesL = new ArrayList<Base>(4);
	public ArrayList<Base> subBasesR = new ArrayList<Base>(4);
	/**
	 * @param baseP 
	 */
	public CopyOfBaseRoot(String basesP) {
		this.bases = basesP;
		quantity++;
		createBases();
	}
	
	/**
	 * @param leftP
	 */
	private void createBases() {
		for (int i = 0; i < Functions.bases.length(); i++) {
			char base = Functions.bases.charAt(i);
			String subString = bases.substring(0, bases.length() - 2) + base;
			subBasesL.add(new Base(subString, true));
		}
		for (int i = 0; i < Functions.bases.length(); i++) {
			char base = Functions.bases.charAt(i);
			String subString = base + bases.substring(2, bases.length());
			subBasesR.add(new Base(subString, false));
		}
	}
	/**
	 * @param base
	 */
	public void addBase(String base) {
		quantity++;
		if (base.length() > 1) {
			char next;
			next = base.charAt(base.length()-2);
			next = base.charAt(2);
			int baseIdx =Functions.mapNumber.get(next);
		}
	}
}
