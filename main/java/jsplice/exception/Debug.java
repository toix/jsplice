/**
 * 
 */
package jsplice.exception;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Debug extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * 
	 */
	public Debug() {
		super();
	}

	/**
	 * @param message
	 */
	public Debug(String message) {
		super(message);
	}

	/**
	 * @param cause
	 */
	public Debug(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message
	 * @param cause
	 */
	public Debug(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * @param message
	 * @param cause
	 * @param enableSuppression
	 * @param writableStackTrace
	 */
	public Debug(String message, Throwable cause,
			boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}
}
