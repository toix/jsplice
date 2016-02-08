/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.Arrays;

//	uint bin;            	"Indexing field to speed chromosome range queries"
//	string name;        	"Name of gene (usually transcript_id from GTF)"
//	string chrom;       	"Chromosome name"
//	char[1] strand;     	"+ or - for strand"
//	uint txStart;       	"Transcription start position"
//	uint txEnd;         	"Transcription end position"
//	uint cdsStart;      	"Coding region start"
//	uint cdsEnd;        	"Coding region end"
//	uint exonCount;     	"Number of exons"
//	uint[exonCount] exonStarts; "Exon start positions"
//	uint[exonCount] exonEnds;   "Exon end positions"
//	uint score;           "score"
//	string name2;       	"Alternate name (e.g. gene_id from GTF)"
//	string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
//	string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
//	lstring exonFrames; 	"Exon frame {0,1,2}, or -1 if no frame for exon"

/**
 * Information of a line from the RefGene file parsed into an Object
 * 
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class RefGene {

	/**
	 * Names of the RefGene file titles
	*/
	private final String[] identifiers = {"bin", "name", "chrom", "strand",
			"txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
			"exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
			"cdsEndStat", "exonFrames" };
	
	private String refSeqAccession;

	private String chromosome;
	/**
	 * Strand Orientation of the mRNA
	 * True - plus strand
	 * False - minus strand
	*/
	private boolean isPlusStrand;

	private int transcriptionStart;

	private int transcriptionEnd;

	private int exonCount;
	/**
	 * Exon start positions on the + strand of the chromosome
	 * Exon end positions on the - strand of the chromosome
	 */
	private ArrayList<Integer> exonStarts = new ArrayList<Integer>();
	/**
	 * Exon end positions on the + strand of the chromosome
	 * Exon start positions on the - strand of the chromosome
	 */
	private ArrayList<Integer> exonEnds = new ArrayList<Integer>();

	/**
	 * GTF ID of the Gene
	 */
	private String gene;

	/**
	 * Constructs an Object from a RefGene line and validates Data
	 * @param line Space separated RefGene line split into an Array
	 * @throws IllegalArgumentException
	 */
	public RefGene(String[] line) throws IllegalArgumentException {
		this(line, true);
	}
	
	/**
	 * Constructs an Object from a RefGene line and validates Data
	 * @param line Space separated RefGene line split into an Array
	 * @param validateChromosome If false all chromosome names will be accepted
	 * @throws IllegalArgumentException
	 */
	public RefGene(String[] line, boolean validateChromosome) throws IllegalArgumentException {
		setRefSeqAccession(line[1]);
		setChromosome(line[2], validateChromosome);
		setPlusStrand(line[3]);
		setTranscriptionStart(line[4]);
		setTranscriptionEnd(line[5]);
		setExonCount(line[8]);
		setExonStarts(line[9]);
		setExonEnds(line[10]);
		setGeneName(line[12]);
	}
	
	/**
	 * Checks whether there is a exon start or exon end next to the junction
	 * @param variantPosition
	 * @return
	 */
	public boolean hasJunctionNextTo(int variantPosition){
		boolean result = false;
		int range = Config.getLengthTrainingSequence();
		for (Integer geneJunction : exonStarts) {
			if (geneJunction+range > variantPosition && geneJunction - range < variantPosition) {
				result = true;
			}
		}
		for (Integer geneJunction : exonEnds) {
			if (geneJunction+range > variantPosition && geneJunction - range < variantPosition) {
				result = true;
			}
		}
		return result;
	}
	
	
	private void setRefSeqAccession(String accessionParam){
		String accession = accessionParam;
		if(accession.contains("."))
			accession = accession.split("\\.", 2)[0];
		Integer.parseInt(accession.substring(3));
		if(!accession.subSequence(0, 3).equals("NM_"))
			throw new IllegalArgumentException("Illegal RefSeq accession number:\t" + accessionParam);
		this.refSeqAccession = accession;
	}
	private void setChromosome(String chromosome, boolean validate) {
		if (chromosome.length() > 2 && chromosome.substring(0, 3).equals("chr")) {
			chromosome = chromosome.substring(3, chromosome.length());
		}
		if (validate && !Arrays.asList(Config.getChromosomenames()).contains(chromosome)) {
			throw new IllegalArgumentException("Illegal chromosome:\t" + chromosome);
		}
		this.chromosome = chromosome;
	}
	private void setPlusStrand(String strand) throws IllegalArgumentException{
		if (strand.equals("+")){
			isPlusStrand = true;
		} else if (strand.equals("-")) {
			isPlusStrand = false;
		} else {
			throw new IllegalArgumentException("Invalid Strand orientation:\t" + strand);
		}
	}
	private void setTranscriptionStart(int txStart){
		this.transcriptionStart = txStart;
	}
	private void setTranscriptionStart(String txStart){
		setTranscriptionStart(Integer.parseInt(txStart));
	}
	
	private void setTranscriptionEnd(int txEnd){
		this.transcriptionEnd = txEnd;
	}
	private void setTranscriptionEnd(String txEnd){
		setTranscriptionEnd(Integer.parseInt(txEnd));
	}
	private void setExonCount(String exonCount){
		this.exonCount = Integer.parseInt(exonCount);
	}

	private void setExonStarts(String exonStartStrings){
		String[] exonStartSplit = exonStartStrings.split(",");
		for (int i = 0; i < exonStartSplit.length; i++) {
			exonStarts.add(Integer.parseInt(exonStartSplit[i]));
		}
	}
	private void setExonEnds(String exonEndStrings){
		String[] exonEndSplit = exonEndStrings.split(",");
		for (int i = 0; i < exonEndSplit.length; i++) {
			exonEnds.add(Integer.parseInt(exonEndSplit[i]));
		}
	}


	/**
	 * @param string
	 */
	private void setGeneName(String geneName) {
		this.gene = geneName;
		
	}

	public String[] getIdentifiers() {
		return identifiers;
	}


	public String getRefSeqAccession() {
		return refSeqAccession;
	}


	public String getChromosome() {
		return chromosome;
	}


	public boolean isPlusStrand() {
		return isPlusStrand;
	}


	public int getTranscriptionStart() {
		return transcriptionStart;
	}


	public int getTranscriptionEnd() {
		return transcriptionEnd;
	}


	public int getExonCount() {
		return exonCount;
	}


	public ArrayList<Integer> getExonStarts() {
		return exonStarts;
	}


	public ArrayList<Integer> getExonEnds() {
		return exonEnds;
	}

	/**
	 * Add 1 to all exon starts
	 * @param addition
	 */
	public void addOneToExonStarts() {
		for (int i = 0; i < exonStarts.size(); i++) {
			exonStarts.set(i, exonStarts.get(i)+1);
		}
	}
	
	/**
	 * Subtract 1 from all exon ends
	 * @param addition
	 */
	public void subOneFromExonEnds() {
		for (int i = 0; i < exonEnds.size(); i++) {
			exonEnds.set(i, exonEnds.get(i)-1);
		}
	}

	public String getGene() {
		return gene;
	}
}
