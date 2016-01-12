/**
 * This utility class contains convenience methods to detect and use
 * respective RDKit functionality of operating systems.
 * 
 */
public final class SystemUtils {

  //
  // Static Methods
  //

  /**
   * Returns true, if the currently used operating system is Microsoft Windows.
   * 
   * @return True, if Windows. False otherwise.
   */
  public static boolean isWindows() {
    final String os = System.getProperty("os.name").toLowerCase();
    // windows
    return (os.indexOf("win") >= 0);
  }

  /**
   * Returns true, if the currently used operating system is Mac OS.
   * 
   * @return True, if Mac. False otherwise.
   */
  public static boolean isMac() {
    final String os = System.getProperty("os.name").toLowerCase();
    // Mac
    return (os.indexOf("mac") >= 0);
  }

  /**
   * Returns true, if the currently used operating system is Unix or Linux.
   * 
   * @return True, if Unix or Linux. False otherwise.
   */
  public static boolean isUnix() {
    final String os = System.getProperty("os.name").toLowerCase();
    // linux or unix
    return (os.indexOf("nix") >= 0 || os.indexOf("nux") >= 0);
  }

  //
  // Constructor
  //

  /**
   * This constructor serves only the purpose to avoid instantiation of this class.
   */
  private SystemUtils() {
    // To avoid instantiation of this class.
  }
}
