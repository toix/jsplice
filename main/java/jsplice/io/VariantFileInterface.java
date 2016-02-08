/**
 * 
 */
package jsplice.io;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public interface VariantFileInterface {
	
	public void readFile(String variantFileName);
	public void addVariant(String line);
}
