package old;

import java.util.ArrayList;

import com.google.common.primitives.Doubles;

import jsplice.data.Config;
import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.tools.Functions;

/**
 * Gene Object storing Variants
 * @author Tobias Gresser (gresserT@gmail.com)
 * TODO Gene loeschen: Gene dürfen mehrmals vorkommen, nur die gleiche Junction mehrmals sollte vermieden werden
 */
public class Gene {

	/**
	 * The Gene ID
	 */
	private String geneId;
	/**
	 * Variants belonging to the gene
	 */
	private ArrayList<Variant> variants = new ArrayList<Variant>();
	/**
	 * Consensus sequence over all variant sequences
	 */
	private String consensus;
	
	/**
	 * 
	 */
	public Gene(String geneId, Variant variant) {
		this.geneId = geneId;
		variants.add(variant);
		consensus = variant.getSequence().toString();
	}
	
	public void add(Variant variant){
		variants.add(variant);
		calculateConsensus();
	}
	
	private void calculateConsensus(){
		char[] consensus = new char[Config.getLengthModelSequence()];
		// get sequences from variants
		ArrayList<Sequence> sequences = getSequences();
		// calculate probability
		try {
			double[][] probability = Functions.getFrequencies(sequences, true);
			for (int i = 0; i < probability.length; i++) {
				// index of the max
				double max = Doubles.max(probability[i]);
				int maxIdx = Doubles.indexOf(probability[i], max);
				consensus[i] = Functions.bases.charAt(maxIdx);
			}
		} catch (NullPointerException e) {
			Functions.arrayToString((String[]) sequences.toArray());
			e.printStackTrace();
		}
	}
	
	/**
	 * @return
	 */
	private ArrayList<Sequence> getSequences() {
		return null;
	}

	public int size(){
		return variants.size();
	}
	
	public Variant get(int i){
		return variants.get(i);
	}

	public String getGeneId() {
		return geneId;
	}

	public String getConsensus() {
		return consensus;
	}
}
