package old;

import java.io.File;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class FastaReferenceChromosomeFile extends IndexedFastaSequenceFile {
	
	public FastaReferenceChromosomeFile(File faFile, File faiFile){
        super(faFile, new FastaSequenceIndex(faiFile));
	}
}