/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.Collections;

import jsplice.exception.Log;
import jsplice.io.Variants;
import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class VariantsTupel {

	public static final double percentageTest = 0.01;
	public Variants test;
	public Variants train;
	
	/**
	 * 
	 */
	public VariantsTupel() {
		test = new Variants();
		train = new Variants();
	}
	
	/**
	 * @param variantsP
	 * @return Crate two Variants; One with train variants and the other with test variants.
	 */
	public static VariantsTupel createTestPool(Variants variantsP) {
		VariantsTupel variantsTupel = new VariantsTupel();
		// create train and test variants
		int testSize = (int) Functions.round(variantsP.size() * percentageTest, 0);
		// create array with unique random indices
		ArrayList<Integer> idxList = new ArrayList<Integer>();
		for (int i = 0; i < variantsP.size(); i++) {
			idxList.add(new Integer(i));
		}
	    Collections.shuffle(idxList);
	    // separate test and train variants
	    for (int i = 0; i < testSize ; i++) {
	    	variantsTupel.test.add(variantsP.get(idxList.get(i)));
		}
	    for (int i = testSize; i < variantsP.size() ; i++) {
	    	variantsTupel.train.add(variantsP.get(idxList.get(i)));
	    }
		Log.add(variantsP.size() + " variants = " + variantsTupel.train.size() + " train variants + " + variantsTupel.test.size() + " test variants", 3);
		return variantsTupel;
	}

}
