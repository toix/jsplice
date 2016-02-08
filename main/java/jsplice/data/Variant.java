package jsplice.data;

import java.util.Arrays;

import jsplice.exception.Log;
import jsplice.exception.RefGeneNotFoundException;
import jsplice.tools.Functions;

/**
 * Information of a variant line from the ClinVar file parsed into an Object
 * 
 * @author Tobias Gresser (gresserT@gmail.com)
 * TODO change Variant dependent fields in field Objects (sequence) when Variant field is changed
 */
public class Variant {

	/**
	 * list of Gene IDs overlapping the variation
	 */
	String gene;
	/**
	 * name of the reference genome assembly version on which genome locations
	 * are based
	 */
	String assembly;
	/**
	 * number or name of the chromosome
	 */
	String chromosome;
	/**
	 * start location on the chromosome (pter->qter)
	 */
	int start;
	/**
	 * end location on the chromosome (pter->qter)
	 */
	int stop;
	/**
	 * chromosome number and band on the chromosome (ISCN Ideogram)
	 */
	String cytogenetic;
	/**
	 * mRNA-based HGVS expression with RefSeq as reference
	 */
	String dnaHgvs;
	/**
	 * RefSeq accession number
	 */
	String refSeqAccession;
	/**
	 * Version of RefSeq accession number for the mRNA
	 */
	int refSeqVersion;
	/**
	 * Sequence of the reference
	*/
	String ref;
	/**
	 * Sequence of the variant
	*/
	String alt;
	/**
	 * Position of the exon junction on the chromosome
	 */
	int distanceToExon;
	/**
	 * Object with additional information about the RefSeq of corresponding to the variant.
	 */
	RefGene refSeq;
	/**
	 * The reference sequence around the junction next to the variant
	 */
	Sequence sequence;
	/**
	 * Is true if variant is related to an exon start junction
	 */
	boolean acceptor;
	/**
	 * Junction absolute position (first or last base of the exon) next to the variant on the chromosome
	 */
	int junctionPosition;
	/**
	 * index of reference SNP in dbSNP, -1 if not available
	 */
	String snp;
	/**
	 * Is true if the change is expected to create a cryptic splice site
	 */
	boolean cryptic;
	/**
	 * DON'T USE
	 * The constructor proves validity for all relevant Annotations
	*/
	public Variant(String[] annotationLine, boolean useGrch38){}
	
	/**
	 * copy Constructor
	 * @param variant
	 */
	public Variant(Variant variant){
		gene = variant.getGene();
		assembly = variant.getAssembly();
		chromosome = variant.getChromosome();
		start = variant.getStart();
		stop = variant.getStart();
		cytogenetic = variant.getCytogenetic();
		dnaHgvs = variant.getDnaHgvs();
		refSeqAccession = variant.getRefSeqAccession();
		refSeqVersion = variant.getRefSeqVersion();
		ref = variant.getRef();
		alt = variant.getAlt();
		distanceToExon = variant.getDistanceToExon();
		refSeq = variant.getRefSeq();
		sequence = variant.getSequence();
		acceptor = variant.isAcceptorSite();
		junctionPosition = variant.getJunctionPosition();
		snp = variant.getSnp();
		cryptic = variant.isCryptic();
	}
	
	/**
	 * @return
	 */
	public String getCytogenetic() {
		return cytogenetic;
	}
	/**
	 * @param string
	 */
	public void setChromosome(String chromosome) {
		if (chromosome.length() > 2 && chromosome.substring(0, 3).equals("chr")) {
			chromosome = chromosome.substring(3, chromosome.length());
		}
		if (!Arrays.asList(Config.getChromosomenames()).contains(chromosome)) {
			throw new IllegalArgumentException("Illegal chromosome:\t" + chromosome);
		}
		this.chromosome = chromosome;
	}

	/**
	 * Parse the RefSeq mRNA accession number without "NM_" prefix and version
	*/
	public void setRefSeqAccession(String refSeq) {
		if(refSeq.contains("."))
			refSeq = refSeq.split("\\.", 2)[0];
		Integer.parseInt(refSeq.substring(3));
		if(!refSeq.subSequence(0, 3).equals("NM_"))
			throw new IllegalArgumentException("Illegal RefSeq accession number:\t" + refSeq);
		this.refSeqAccession = refSeq;
	}
	
	public void setRefSeqVersion(int refSeqVersion) {
		this.refSeqVersion = refSeqVersion;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}
	
	public void setAlt(String alt) {
		this.alt = alt;
	}

	/**
	 * @param cryptic
	 */
	public void setCryptic(boolean cryptic) {
		this.cryptic = cryptic;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}
	/**
	 * Distance to the next exon end or start <br/>
	 * Is negative if it refers to an exon start
	 * @return
	 */
	public int getDistanceToExon() {
		return distanceToExon;
	}

	public String getDnaHgvs() {
		return dnaHgvs;
	}

	/**
	 * Returns the RefSeq accession number
	*/
	public String getRefSeqAccession() {
		return refSeqAccession;
	}

	/**
	 * @return
	 */
	public int getRefSeqVersion() {
		return refSeqVersion;
	}

	public String getAssembly() {
		return assembly;
	}

	public String getChromosome() {
		return chromosome;
	}
	
	public RefGene getRefSeq() {
		return refSeq;
	}

	public void setRefGene(RefGene refGene) throws RefGeneNotFoundException {
		if (refGene.isPlusStrand()) {
			junctionPosition = start - distanceToExon;
		} else {
			junctionPosition = start + distanceToExon;
		}
		// Exon start positions
		if(acceptor == refGene.isPlusStrand()){
			if (refGene.getExonStarts().contains(junctionPosition)) {
				this.refSeq = refGene;
			} else if (refGene.getExonStarts().contains(junctionPosition-1)) {
				// Addition to exonJunction if RefGene contains junction positions of the intron (instead of the exon)
				Log.add("Exon start positions of RefSeq Gene " + refGene.getRefSeqAccession() + " lay in introns and were corrected.", 1);
				refGene.addOneToExonStarts();
				this.refSeq = refGene;
			} else {
				throw new RefGeneNotFoundException("Variant with cHGVS "
						+ dnaHgvs + " and estimated junction at " + junctionPosition
						+ " doesn't fit to any junction of the RefSeqGene file");
			}
		}
		// Exon end positions
		else {
			if (refGene.getExonEnds().contains(junctionPosition)) {
				this.refSeq = refGene;
			} else if (refGene.getExonEnds().contains(junctionPosition+1)) {
				// Subtraction to exonJunction if RefGene contains junction positions of the intron (instead of the exon)
				Log.add("RefSeqGene " + refGene.getRefSeqAccession() + " exon end positions lay in introns.", 4);
				refGene.subOneFromExonEnds();
				this.refSeq = refGene;
			} else {
				throw new RefGeneNotFoundException("Variant with cHGVS "
						+ dnaHgvs + " and estimated junction at " + junctionPosition
						+ " doesn't fit to any junction of the RefSeqGene file");
			}
		}
	}
	/**
	 * @return The sequence around the junction next to the variant
	 */
	public Sequence getSequence() {
		return sequence;
	}
	
	/**
	 * Set the sequence for the variant
	 * @param sequenceStr The sequence around the junction next to the variant
	 * @param extendedExonLength
	 * @return
	 */
	public Sequence setSequence(String sequenceStr, int extendedExonLength, int length, int exonLength) throws RuntimeException {
		this.sequence = new Sequence(sequenceStr, extendedExonLength, length, exonLength, this);
		return this.sequence;
	}
	
	/**
	 * @return Junction position (first or last base of the exon) next to the variant
	 */
	public int getJunctionPosition() {
		return junctionPosition;
	}

	/**
	 * @return Is true if variant is related to an exon start junction
	 */
	public boolean isAcceptorSite() {
		return acceptor;
	}

	public String getRef() {
		return ref;
	}

	/**
	 * @return Sequence of the variant
	*/
	public String getAlt() {
		return alt;
	}

	public String getLine() {
		return chromosome + '\t' + start + '\t' + refSeqAccession + "."
				+ refSeqVersion + '\t' + ref + '\t' + alt + '\t' + ".";
	}

	public String getGene() {
		return gene;
	}
	
	/**
	 * @return
	 */
	private String getSnp() {
		return snp;
	}

	/**
	 * @return
	 */
	public boolean isCryptic() {
		return cryptic;
	}

	@Override
	public String toString(){
		return "gene: " + gene + " chr: " + chromosome + " start: " + start + " HGVS: " + dnaHgvs + " RefSeq: " + refSeqAccession + " distExon: " + distanceToExon + " inExJ: " + acceptor + " junPos: " + junctionPosition  + "\n\tsequence:" + sequence;
	}

	/**
	 * 
	 */
	public void changeRandomBase() {
		distanceToExon = Functions.random(3, sequence.getLengthIntron());
		if (sequence.isAcceptor()) {
			distanceToExon = -distanceToExon; 
		}
		int positionChangeSequence = sequence.getPositionJunction() + distanceToExon;
		
		alt = Functions.bases.charAt(Functions.random(0, 3)) + "";
		while (alt.charAt(0) == sequence.charAt(positionChangeSequence)) {
			alt = Functions.bases.charAt(Functions.random(0, 3)) + "";
		}
		
//		sequenceString = sequenceString.substring(0, positionChangeAbsolute) + alt + sequenceString.substring(positionChangeAbsolute + 1, length());
		
		setSequence(sequence.getStringExtended(), sequence.getLengthExonExtended(), sequence.lengthExtended() / 2,
				sequence.getLengthExonExtended() / 2);
	}
}
