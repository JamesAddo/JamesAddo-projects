/*MedChemCompanion MMPMaker's Molecule Fragmenter and Matched Molecular Pairs Identifier based on Java implementation of Jameed Hussain and Ceara Rea algorithm-Computationally Efficient Algorithm to Identify Matched Molecular Pairs (MMPs) in Large Data Sets, J. Chem. Inf. Model., 2010, 50 (3), pp 339–348.
/
/*
 * Copyright (C)2015, James Addo.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 * 
 * - Neither the name of James Addo
 *   nor the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.RDKit.GenericRDKitException;
import org.RDKit.Int_Vect;
import org.RDKit.ROMol;

public class RDKit {

  //
  // Constants
  //

  /** The logger instance. */
  private static Logger LOGGER = Logger.getLogger(RDKit.class.getName());

  private static final String OS_WIN32 = "win32";

  private static final String OS_LINUX = "linux";

  private static final String OS_MACOSX = "macosx";

  private static final String ARCH_X86 = "x86";

  private static final String ARCH_X86_64 = "x86_64";

  /** List of libraries to be loaded for different operating systems (lib order is important). */
  private static final Map<String, String[]> LIBRARIES = new HashMap<String, String[]>();

  /** We define here what libraries are necessary to run the RDKit for the different supported platforms. */
  static {
    LIBRARIES.put(OS_WIN32 + "." + ARCH_X86, new String[] { "boost_system-vc100-mt-1_51.dll", "GraphMolWrap.dll" });
    LIBRARIES.put(OS_WIN32 + "." +  ARCH_X86_64, new String[] { "boost_system-vc100-mt-1_51.dll", "GraphMolWrap.dll" });
    LIBRARIES.put(OS_LINUX + "." +  ARCH_X86, new String[] { "libGraphMolWrap.so" });
    LIBRARIES.put(OS_LINUX + "." +  ARCH_X86_64, new String[] { "libGraphMolWrap.so" });
    LIBRARIES.put(OS_MACOSX + "." +  ARCH_X86_64, new String[] { "libGraphMolWrap.jnilib" });
  }


  /** Flag to determine, if an activation was successful. */
  private static boolean g_bActivated = false;

  /** Flag to determine, if an activation was already performed. */
  private static boolean g_bActivationRan = false;

  //
  // Constructor
  //

  private RDKit() {
    // Only here to avoid instantiation of this utility class
  }

  //
  // Static Public Methods
  //

  /**
   * Activates the RDKit to be used from a Java Application. This activation
   * call gets only executed once. If it fails there is no recovery possible
   * for the current Java VM. Subsequent calls to this method result in the
   * same boolean result as the first call. The native libraries are expected
   * to be found in lib/native/os/<OS Dependent Path>. This path gets resolved
   * to an absolute path by the Java VM (File class).
   * 
   * @return True, if activation was successful. False otherwise.
   */
  public static synchronized boolean activate() {
    return activate(null);
  }

  /**
   * Activates the RDKit to be used from a Java Application. This activation
   * call gets only executed once. If it fails there is no recovery possible
   * for the current Java VM. Subsequent calls to this method result in the
   * same boolean result as the first call.
   * 
   * @param strPath
   *            An absolute or relative path to a directory, which contains OS
   *            dependent sub directories win32, macosx and linux, which again
   *            contain OS architecture dependent sub directories x86 and
   *            x86_64. If the path is relative it gets resolved to an
   *            absolute path by the Java VM (File class).
   * 
   * 
   * @return True, if activation was successful. False otherwise.
   */
  public static synchronized boolean activate(final String strPath) {
    if (!g_bActivationRan) {
      g_bActivationRan = true;

      try {
        String strRelativePath = (strPath == null ? "lib/native/os/"
            : strPath);

        // Determine operating system
        String strOS;
        if (SystemUtils.isWindows()) {
          strOS = OS_WIN32;
        }
        else if (SystemUtils.isMac()) {
          strOS = OS_MACOSX;
        }
        else if (SystemUtils.isUnix()) {
          strOS = OS_LINUX;
        }
        else {
          throw new UnsatisfiedLinkError(
              "Operating system is not supported from RDKit Native Libraries.");
        }

        // Determine OS architecture (32 or 64 bit)
        String strArch;
        if (System.getProperty("os.arch").indexOf("64") != -1) {
          strArch = ARCH_X86_64;

        }
        else {
          strArch = ARCH_X86;
        }

        strRelativePath += strOS + "/" + strArch + "/";
        final String[] arrLibraries = LIBRARIES.get(strOS + "." + strArch);

        if (arrLibraries == null) {
          throw new UnsatisfiedLinkError("Unsupported operating system or architecture.");
        }
        else {
          // Load libraries
          for (final String strLibName : arrLibraries) {
            System.load(new File(strRelativePath, strLibName).getAbsolutePath());
          }
        }

        g_bActivated = true;
      }
      catch (final UnsatisfiedLinkError e) {
        LOGGER.log(Level.SEVERE, "Unable to load RDKit Native Libraries.", e);
      }
    }

    return g_bActivated;
  }
