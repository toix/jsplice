/**
 * 
 */
package jsplice.data;

import static org.junit.Assert.*;

import java.io.IOException;

import jsplice.exception.Log;
import jsplice.io.ClinVarFile;
import jsplice.io.FastaReferenceFiles;
import jsplice.io.RefSeqGeneFile;
import jsplice.io.Variants;
import jsplice.tools.Functions;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class SequenceTest {

	static Sequence sequence;
	/**
	 * @throws IOException 
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws IOException {
		new Config("testData/variant_summary.txt", "", "", "testData/chromFa/chr*.fa", "testdata/RefSeqGenesGRCh37.tsv", "test_log/sequence.log", 25, 7 , false);
		new Functions();
		new Log();
		RefSeqGeneFile refGeneFile = new RefSeqGeneFile(Config.getRefGeneFileName());
		FastaReferenceFiles faRefFile = new FastaReferenceFiles(Config.getChrFileNames());
		ClinVarFile clinVarFilePatho = new ClinVarFile(Config.getClinVarFileName(), false);
		clinVarFilePatho.addRefSeqData(refGeneFile, faRefFile);
		ClinVarFile clinVarFileBen = new ClinVarFile(Config.getClinVarFileName(), true);
		clinVarFileBen.addRefSeqData(refGeneFile, faRefFile);
		Variants clinVarFile = ClinVarFile.concat(clinVarFileBen, clinVarFilePatho);
		Log.setFileLogLevel(1);
		Log.writeToFile();
		Log.setPrintLogLevel(4);
//		System.out.println(Log.toStringStatic());
		sequence = clinVarFile.get(0).getSequence();
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
		
	}

	@Test
	public void testGetJunctionPosition() {
		boolean result = sequence.getPositionJunction() == 18;
//		System.out.println("jun: " + sequence.getJunctionPosition());
		assertTrue(result);
	}
	
	@Test
	public void testGetVariantPosition() {
		boolean result = sequence.getPositionChange() == 3;
//		System.out.println("change: " + sequence.getChangePosition());
		assertTrue(result);
	}
	
	@Test
	public void testGetReference() {
		boolean result = sequence.getReference().equals("T");
		assertTrue(result);
	}
	
	/**
	 * variant pos 43248503
	 * junction pos 43248488
	 * 43248488/50 = 864969 R 38 
	 * -> line 864971 row 23 
	 * TCAGAGGTTCTGCC CTGCTCGTTGGTTCAGAGAAGCAAAAAGATCAGGCA
	 * TCAGAGGTTCTGCCCTGCTCGTTGGTTCAGAGAAGCAAAAAGATCAGGCA
	 * TGCCTGATCTTTTTGCTTCTCTGAACCAACGAGCAGGGCAGAACCTCTGA
	 * TGCCTGATCTTTTTGCTTCTCTGAACCAACGAGCAG GGCAGAACCTCTGA
	 */
	@Test
	public void testToString() {
		boolean result = sequence.toString().equals("CTCTGAACCAACGAGCAG_GGCAGAA junPos: 18 varPos: -15 alt: A intronExon: true");
//		System.out.println(sequence.toString());
		assertTrue(result);
	}

	@Test
	public void testGetAlternateSequence() {
		boolean result = sequence.getStringAlternate().equals("CTCAGAACCAACGAGCAGGGCAGAA");
//		System.out.println(sequence.getAlternateSequence());
		assertTrue(result);
	}
	
	@Test
	public void testGetIntronicPart() {
		boolean result = sequence.getStringIntron().equals("CTCTGAACCAACGAGCAG");
//		System.out.println(sequence.getIntronicPart());
		assertTrue(result);
	}
	
	@Test
	public void testCharAt() {
		boolean result = sequence.charAt(12) == 'G' && sequence.charAt(2) == 'C';
//		System.out.println(sequence.charAt(12) + " " + sequence.charAt(2));
		assertTrue(result);
	}
	
	@Test
	public void testSubstring() {
		boolean result = sequence.substring(2, 5).equals("CTG");
//		System.out.println(sequence.substring(-1, 1));
		assertTrue(result);
	}
	
	@Test
	public void testSubstringRef() {
		boolean result = sequence.substring(-1, 1, false).equals("TC");
//		System.out.println(sequence.substring(-1, 1));
		assertTrue(result);
	}
	
	@Test
	public void testSubstringAlt() {
		boolean result = sequence.substring(2, 5, false).equals("CAG");
//		System.out.println(sequence.substring(2, 5, false));
		assertTrue(result);
	}
}
