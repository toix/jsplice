/**
 * 
 */
package jsplice.tools;

import static org.junit.Assert.*;
import jsplice.data.Config;
import jsplice.exception.Log;
import jsplice.io.ClinVarFile;
import jsplice.io.FastaReferenceFiles;
import jsplice.io.HgmdFile;
import jsplice.io.RefSeqGeneFile;
import jsplice.io.Variants;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class FunctionsTest {

	private static Variants variants;
	private static Sequences sequences;
	private static boolean intronExon;
	
	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		new Config("testData/variant_summary.txt", "testData/hgmd.txt", "", "testData/chromFa/chr*.fa", "testdata/RefSeqGenesGRCh37.tsv", "test_log/sequences.log", 25, 7 , false);
		new Functions();
		new Log();
		intronExon = true;
		
		ClinVarFile clinVariants = new ClinVarFile(Config.getClinVarFileName(), false);
		HgmdFile hgmdVariants = new HgmdFile(Config.getHgmdFileName());
		RefSeqGeneFile refGeneFile = new RefSeqGeneFile(Config.getRefGeneFileName());
		FastaReferenceFiles faRefFile = new FastaReferenceFiles(Config.getChrFileNames());
		clinVariants.addRefSeqData(refGeneFile, faRefFile);
		hgmdVariants.addRefSeqData(refGeneFile, faRefFile);
		Variants variants = Variants.concat(clinVariants, hgmdVariants);
		Variants variantsBeforeExon = Variants.filterVariants(variants, intronExon);
		sequences = new Sequences(variantsBeforeExon, true);
//		System.out.println(sequences);
		
		Log.setFileLogLevel(1);
		Log.writeToFile();
	}

	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		Log.writeToFile();
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void validComplement(){
		String result = Functions.complement("aAcCgtACtGtTNcgNc");
		assertTrue(result.equals("tTgGcaTGaCaANgcNg"));
	}
	
	@Test(expected = IllegalArgumentException.class)
	public void invalidComplement() {
		Functions.complement("acgtACrGTN");
	}
	
	@Test
	public void validReverseComplement() {
		String result = Functions.reverseAndComplement("aAcCgtACtGtTNcgNc");
		assertTrue(result.equals("gNcgNAaCaGTacGgTt"));
	}

	@Test(expected = IllegalArgumentException.class)
	public void invalidReverseComplement() {
		Functions.reverseAndComplement("acgtACäGTN");
	}
	
	@Test
	public void Reverse(){
		String result = Functions.reverse("aAcCgtACtGtTNcgNc");
		assertTrue(result.equals("cNgcNTtGtCAtgCcAa"));
	}
	
	@Test(expected = IllegalArgumentException.class)
	public void invalidReverse() {
		Functions.reverse("acgtACrGTN");
	}
	
	@Test
	public void transposeArray(){
		String[] strArray = {"abcdefg", "bcdefgh", "cdefghi", "defghij"};
		String[] transStrArray = {"abcd", "bcde", "cdef", "defg", "efgh", "fghi", "ghij"};
		String[] result = Functions.transpose(strArray);
		for (int i = 0; i < transStrArray.length; i++) {
			assertTrue(transStrArray[i].equals(result[i]));
		}
	}
	
	@Test
	public void transposeSequence(){
		
		String[] expectation = Functions.transpose(new String[]{sequences.get(0).getStringReference(), sequences.get(1).getStringReference(), sequences.get(2).getStringReference(), sequences.get(3).getStringReference(), sequences.get(4).getStringReference()});
		String[] result = Functions.transpose(sequences.getSequences(), true);
//		System.out.println(Functions.arrayToString(expectation));
//		System.out.println(Functions.arrayToString(result));
		for (int i = 0; i < expectation.length; i++) {
			assertTrue(result[i].equals(expectation[i]));
		}
	}

	/**
	 * Test method for {@link jsplice.tools.Sequences#convertToUpperBases(java.lang.String)}.
	 */
	@Test
	public void testConvertToUpperBases() {
		String result = Functions.convertToUpperBases("AcGttAcgtGGCCA");
		assertTrue(result.equals("ACGTTACGTGGCCA"));
	}

	/**
	 * Test method for {@link jsplice.tools.Sequences#convertToUpperBases(java.lang.String)}.
	 */
	
	@Test (expected = IllegalArgumentException.class)
	public void testBadConvertToUpperBases() {
		Functions.convertToUpperBases("AcGttAcgtGNCCA");
	}
	
	/**
	 * Test method for {@link jsplice.tools.Sequences#convertToUpperBases(java.lang.String)}.
	 */
	
	@Test
	public void testSum() {
		double result = Functions.sum(new double[]{1.2, 2.2, 3.2, 1.1, 5.1, 3.3});
		assertTrue(result == 16.1);
	}

	/**
	 * Test method for {@link jsplice.tools.Sequences#getFrequencies(java.lang.String)}.
	 */
	@Test
	public void testGetProbabilities() {
		double[] result = Functions.getFrequencies("AACCGAGTAC");
		assertArrayEquals(new double[]{0.4, 0.3, 0.2, 0.1}, result, 0);
	}
	
	/**
	 * Test method for {@link jsplice.tools.Sequences#getFrequencies(java.lang.String)}.
	 */
	@Test
	public void testGetSequenceProbabilities() {
		double[][] probabilities = Functions.getFrequencies(sequences.getSequences(), true);
		String[] transposed = Functions.transpose(sequences.getSequences(), true);
		
		for (int i = 0; i < transposed.length; i++) {
			double[] expectation = Functions.getFrequencies(transposed[i]);
			boolean result = Functions.equals(probabilities[i], expectation);
//			System.out.println(transposed[i]);
//			System.out.println(Functions.arrayToString(expectation));
//			System.out.println(Functions.arrayToString(probabilities[i]));
			assertTrue(result);
		}
	}
	

	/**
	 * Test method for {@link jsplice.tools.Sequences#getFrequencies(java.lang.String)}.
	 */
	@Test (expected = IllegalArgumentException.class)
	public void testBadGetProbabilities() {
		Functions.getFrequencies("AACCNTGTAC");
	}
}