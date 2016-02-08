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
public class Base {

	public String bases;
	public int quantity = 0;
	public boolean left;
	public ArrayList<Base> subBases = new ArrayList<Base>(4);
	/**
	 * @param baseP 
	 * 
	 */
	public Base(String basesP, boolean leftP) {
		this.bases = basesP;
		quantity++;
		this.left = leftP;
		createBases();
	}
	/**
	 * @param left
	 */
	private void createBases() {
		if (left) {
			for (int i = 0; i < Functions.bases.length(); i++) {
				char base = Functions.bases.charAt(i);
				String subString = bases.substring(0, bases.length()-2) + base;
				subBases.add(new Base(subString, left));
			}
		} else {
			for (int i = 0; i < Functions.bases.length(); i++) {
				char base = Functions.bases.charAt(i);
				String subString = base + bases.substring(2, bases.length());
				subBases.add(new Base(subString, left));
			}
		}
		
	}
	/**
	 * @param base2
	 */
	public void addBase(String base) {
		quantity++;
		if (base.length() > 1) {
			char next;
			if (left) {
				next = base.charAt(base.length()-2);
			} else {
				next = base.charAt(2);
			}
			int baseIdx =Functions.mapNumber.get(next);
		}
	}

}
