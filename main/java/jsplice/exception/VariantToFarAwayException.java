/**
 * 
 */
package jsplice.exception;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class VariantToFarAwayException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * 
	 */
	public VariantToFarAwayException() {
		super();
	}

	/**
	 * @param message
	 */
	public VariantToFarAwayException(String message) {
		super(message);
	}

	/**
	 * @param cause
	 */
	public VariantToFarAwayException(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message
	 * @param cause
	 */
	public VariantToFarAwayException(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * @param message
	 * @param cause
	 * @param enableSuppression
	 * @param writableStackTrace
	 */
	public VariantToFarAwayException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}

}
