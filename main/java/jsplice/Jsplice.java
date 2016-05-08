package jsplice;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import jsplice.data.Config;
import jsplice.exception.Log;
import jsplice.io.ClinVarFile;
import jsplice.io.FastaReferenceFiles;
import jsplice.io.HgmdFile;
import jsplice.io.RefSeqGeneFile;
import jsplice.io.Variants;
import jsplice.tools.AlgorithmAdministrator;
import jsplice.tools.Functions;

public class Jsplice {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws Exception 
	 */
	public static void main(String[] args) throws IOException  {
		try {
		new Config();
		new Functions();
		new Log();
		
		
		// RefSeq
		Log.add("\n - - - Reference - - - ", 3);
		FastaReferenceFiles faRefFile = new FastaReferenceFiles(Config.getChrFileNames());
		RefSeqGeneFile refGeneFile = new RefSeqGeneFile(Config.getRefGeneFileName());
		
		int lenIntron = Config.getLengthModelIntron();
		System.gc();
		Config.setLengthModelIntron(lenIntron);
		Log.add("len intron: " + Config.getLengthModelIntron());
		// ClinVar pathogenic
		Log.add("\n - - - ClinVar pathogene - - - ", 3);
		ClinVarFile clinVarFile = new ClinVarFile(Config.getClinVarFileName(), false);
		clinVarFile.addRefSeqData(refGeneFile, faRefFile);
		// HGMD (pathogene)
		Log.add("\n - - - HGMD pathogene - - - ", 3);
		HgmdFile hgmdFile = new HgmdFile(Config.getHgmdFileName());
		hgmdFile.addRefSeqData(refGeneFile, faRefFile);
		// pathogenic
		Variants varFile = Variants.concat(clinVarFile, hgmdFile);
  		// Benign
		Log.add("\n - - - ClinVar benign - - - ", 3);
		ClinVarFile clinVarFileBenign = new ClinVarFile(Config.getClinVarFileName(), true);
		clinVarFileBenign.addRefSeqData(refGeneFile, faRefFile);
		// prediction
		AlgorithmAdministrator algoAdmin = new AlgorithmAdministrator(varFile, clinVarFileBenign);
		
		
		//vcf
		VCFFileReader vcfReader = new VCFFileReader(new File(Config.getVcfFileName()), false);
        CloseableIterator<VariantContext> iterator = vcfReader.iterator();
//		if (iterator.hasNext()) {
//			System.out.println(iterator.next().getChr());
//		}
        vcfReader.close();
        Log.writeToFile();
		Log.close();
		} catch (Exception e) {
		  try {
			Log.writeToFile();
			Log.close();
			e.printStackTrace();
		  } catch (Exception e2) {
		    Log.close();
            e.printStackTrace();
		  }
		}
	}
}