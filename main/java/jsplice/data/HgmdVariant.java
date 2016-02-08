/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import jsplice.exception.DebuggingException;
import jsplice.exception.Log;
import jsplice.exception.VariantToFarAwayException;
import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class HgmdVariant extends Variant {

	/**
	 * Disease or condition associated with the mutation
	 */
	String disease;
	/**
	 * Full name of the gene
	 */
	String geneName;
	/**
	 * Identifier of the GDB Genome Database
	 */
	String gdbId;
	/**
	 * Identifier if the OMIM database
	 */
	String omim;
	/**
	 * Evidence of the association between disease and mutation <br/>
	 * DM - disease-causing mutations <br/>
	 * DP - polymorphism with significant association; disease related function assumed <br/>
	 * DFP - polymorphism with significant association; disease related function <br/>
	 * FP - polymorphism with significant association to gene, no disease association <br/>
	 */
	String evidence;
	/**
	 * PubMed id for the reference
	 */
	String pubMedId;
	/**
	 * Comment of the curator
	 */
	String comment;
	/**
	 * HGMD accession
	 */
	String hgmdAccession;
	/**
	 * Date when the mutation was added to the database
	 */
	String date;
	/**
	 * Type of the change: Deletion, Mutation, Splice, Promotinal, Insertion, X indel, Grosdel, E repeat, N big indel
	 */
	char type;
	/**
	 * Pattern for Substitutions
	 * NM_ - Prefix of RefSeq mRNA accession number
	 */
	final Pattern checkSubstitution = Pattern.compile("[1-9][0-9]*[+-][1-9][0-9]*([ACGT]>[ACGT])$");
	
	/**
	 * @param lines
	 * @param annotationTitles
	 * @param irrelevantAnnotations
	 * @param useGrch38
	 * @throws NumberFormatException
	 * @throws IllegalArgumentException
	 */
	public HgmdVariant(String[] lines) throws NumberFormatException, IllegalArgumentException {
		super(lines, false);
		// 29 hgmd accession
		hgmdAccession = lines[28];
		if (!Pattern.matches("C[DMSPIXGEN][0-9]{6,7}", hgmdAccession)) {
			throw new IllegalArgumentException("Illegal HGMD accession:\t\"" + hgmdAccession + "\"");
		}
		// if a relevant annotation is empty or contains no data, throw IllegalArgumentException
		List<Integer> relevantAnnotations = Arrays.asList(1, 12, 15, 16, 17, 18, 28, 30);
		for (int i = 0; i < lines.length; i++) {
			if (relevantAnnotations.contains(i)	&& Pattern.matches("^[\\s\\-]*$", lines[i])){
				throw new IllegalArgumentException("The following row contains no data: " + i + "\t HGMD:" + hgmdAccession);
			}
		}
		// 12 DNA change
		dnaHgvs = lines[12];
		// check for single nucleotide splice variant or throw exception
		if(!checkSubstitution.matcher(dnaHgvs).matches())// && !checkOtherHGVS.matcher(annotationLine[18]).matches())
			throw new IllegalArgumentException("Illegal DNA change or non-intronic mutation: " + dnaHgvs + "\t HGMD: " + hgmdAccession);
		// 0 disease
		disease = lines[0];
		// 1 gene
		gene = lines[1].split(";")[0];
		// 21 chromosome
		setChromosome(lines[15]);
		// 2 chromosome band
		setCytogenic(lines[2]);
		// 3 gene ID
		geneName = lines[3];
		// 4 gdbId
		gdbId = lines[4];
		// 5 Omim
		omim = lines[5];
		// 20 SNP
		snp = lines[14];
		// 22 start
		start =Integer.parseInt(lines[16]);
		// 23 stop
		stop = Integer.parseInt(lines[17]);
		// 24 evidence
		evidence = lines[18];
//		if(evidence.equals("DM?"))
//			evidence = "DP";
		final ArrayList<String> checkOrigin = new ArrayList<String>(Arrays.asList("DM","DP","DFP"));
		if (!checkOrigin.contains(evidence)) {
			throw new IllegalArgumentException("Illegal evidence: " + evidence + "\t HGMD: " + hgmdAccession);
		}
		// 25 pubMed
		pubMedId =  lines[25];
		// 27 comment
		comment = lines[27];
		// 29 date
		date = lines[29];
		// 30 type
		type = lines[30].charAt(0);
		if(type != 'S')
			throw new IllegalArgumentException("Non splice DNA variant: \"" + type + "\"\t HGMD: " + hgmdAccession);
		// Parse additional Information from File
		parseInformation();
		Log.add("HGMD Variant with accession " + hgmdAccession + " added.", 1);
	}
	
	/**
	 * Validate and set cytogenetic
	 * @param chrBand chromosome number and band on the chromosome
	 */
	private void setCytogenic(String chrBand) {
		if (!Pattern.matches("^"+ chromosome +"[pq][1-4][0-9](\\.[1-9][0-9]?(\\.[1-9])?)?"
				+ "(-(("+ chromosome +")?[pq])?[1-4][0-9](\\.[1-9][0-9]?(\\.[1-9])?)?)?$",chrBand))
			throw new IllegalArgumentException("Illegal chromosome band:\t" + chrBand + "\t HGMD: " + hgmdAccession);
		cytogenetic =chrBand;
	}

	/**
	 * Parse the Information of the Variant in other fields.
	*/
	private void parseInformation(){
		// before exon start?
		acceptor = dnaHgvs.contains("-");
		// distanceToExon
		int begin;
		int end = dnaHgvs.length()-3;
		if(acceptor){
			begin = dnaHgvs.indexOf('-');
		}
		else{
			begin = dnaHgvs.indexOf('+');
		}
		if(begin == -1)
			throw new DebuggingException("Can't find + or - in DNA change.\t HGMD:" + hgmdAccession);
		if(start != stop)
			throw new DebuggingException("start != stop\tHGMD:" + hgmdAccession);
		distanceToExon = Integer.parseInt(dnaHgvs.substring(begin, end));
		if(Math.abs(distanceToExon) > Config.getLengthTrainingIntron())
			throw new VariantToFarAwayException("The variant with HGVS " + dnaHgvs + " is to far away from junction: " + Math.abs(distanceToExon));
		// alt
		if(!dnaHgvs.contains(">")){
			throw new DebuggingException("DNA change contains no \">\"\t HGMD:" + hgmdAccession);
		}
		String chgvs[] = dnaHgvs.split(">");
		alt = chgvs[1];
		if (!Functions.isValidDNA(alt)) {
			throw new DebuggingException("DNA change contains invalid DNA:" + alt + "\t HGMD:" + hgmdAccession);
		}
		// ref
		chgvs = chgvs[0].split("[0-9]");
		ref = chgvs[chgvs.length-1];
		if (!Functions.isValidDNA(ref)) {
			throw new DebuggingException("DNA change contains invalid DNA:" + ref + "\t HGMD:" + hgmdAccession);
		}
	}
}
