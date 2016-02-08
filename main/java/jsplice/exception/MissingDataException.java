/**
 * 
 */
package jsplice.exception;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class MissingDataException extends RuntimeException {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * 
	 */
	public MissingDataException() {
		super();
	}

	/**
	 * @param message
	 */
	public MissingDataException(String message) {
		super(message);
	}

	/**
	 * @param cause
	 */
	public MissingDataException(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message
	 * @param cause
	 */
	public MissingDataException(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * @param message
	 * @param cause
	 * @param enableSuppression
	 * @param writableStackTrace
	 */
	public MissingDataException(String message, Throwable cause,
			boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}
}