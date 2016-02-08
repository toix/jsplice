/**
 * 
 */
package old;

import java.util.ArrayList;
import java.util.HashMap;

import jsplice.data.Variant;
import jsplice.io.Variants;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 * TODO Genes loeschen: Gene dürfen mehrmals vorkommen, nur die gleiche Junction mehrmals sollte vermieden werden 
 */
public class Genes {

	/**
	 * The Genes containing the Variants
	 */
	ArrayList<Gene> genes = new ArrayList<Gene>();
	/**
	 * Hash the genes by Gene ID
	 */
	HashMap<String, Gene> hashByGeneId = new HashMap<String, Gene>();
	
	/**
	 * 
	 * @param variants
	 */
	public Genes(Variants variants) {
		for (int i = 0; i < variants.size(); i++) {
			add(variants.get(i));
		}
	}
	
	public boolean add(Variant variant){
		String geneId = variant.getGene();
		if (!hashByGeneId.containsKey(geneId)) {
			Gene gene = new Gene(geneId, variant);
			genes.add(gene);
			hashByGeneId.put(geneId, gene);
		} else{
			hashByGeneId.get(geneId).add(variant);
		}
		return true;
	}

	public int size(){
		return genes.size();
	}

	public Gene get(int i){
		return genes.get(i);
	}

	public ArrayList<Gene> getGenes() {
		return genes;
	}
}
