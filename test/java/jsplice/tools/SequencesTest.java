/**
 * 
 */
package jsplice.tools;

import static org.junit.Assert.assertTrue;
import jsplice.data.Config;
import jsplice.exception.Log;
import jsplice.io.ClinVarFile;
import jsplice.io.FastaReferenceFiles;
import jsplice.io.HgmdFile;
import jsplice.io.RefSeqGeneFile;
import jsplice.io.Variants;

import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class SequencesTest {
	
	private static Model model;

	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		new Config("testData/variant_summary.txt", "testData/hgmd.txt", "", "testData/chromFa/chr*.fa", "testdata/RefSeqGenesGRCh37.tsv", "test_log/sequences.log", 25, 7 , false);
		new Functions();
		new Log();
		
		ClinVarFile clinVariants = new ClinVarFile(Config.getClinVarFileName(), false);
		HgmdFile hgmdVariants = new HgmdFile(Config.getHgmdFileName());
		RefSeqGeneFile refGeneFile = new RefSeqGeneFile(Config.getRefGeneFileName());
		FastaReferenceFiles faRefFile = new FastaReferenceFiles(Config.getChrFileNames());
		clinVariants.addRefSeqData(refGeneFile, faRefFile);
		hgmdVariants.addRefSeqData(refGeneFile, faRefFile);
		Variants variants = Variants.concat(clinVariants, hgmdVariants);
		Variants variantsBeforeExon = Variants.filterVariants(variants, true);
		model = new Model(variantsBeforeExon, true);
//		System.out.println(sequences);
		
		Log.setFileLogLevel(1);
		Log.writeToFile();
	}
	
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		Log.writeToFile();
	}
	
	@Before
	public void setUp() throws Exception {
	}
	
	@Test
	public void testGetCorrelation(){
		double correlation = model.getCorrelation(30, 33);
		boolean result = correlation == 0.2132993161855452;
		assertTrue(result);
	}
}
