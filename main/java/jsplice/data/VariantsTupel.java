/**
 * 
 */
package jsplice.data;

import jsplice.io.Variants;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class VariantsTupel {

	public static final double percentageTest = 0.05;
	public Variants test;
	public Variants train;
	
	public VariantsTupel() {
	  super();
	  test = new Variants();
	}

	/**
	 * @param variantsP
	 * @return Crate two Variants; One with train variants and the other with test variants.
	 */
	public static VariantsTupel createTestPool(Variants variantsP, int v) {
		VariantsTupel variantsTupel = new VariantsTupel();
		variantsTupel.train = new Variants(variantsP);
		variantsTupel.test.add(variantsTupel.train.get(v));
		variantsTupel.train.remove(v);
		return variantsTupel;
	}

}
