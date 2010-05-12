/**
 * Fragmentor.x10
 *
 * Basic implementation of MTA in X10, for large scale HF calculation. 
 * Lots of code derived from MeTA Studio and MTA-GAMESS written by me ;)
 * This is a straight implementation, so no plugins are provided as 
 * with MeTA Studio.
 *
 * Primary reference:
 *   [MTA06] V. Ganesh, R. K. Dongare, P. Balanarayan, and S. R. Gadre, 
 *           J. Chem. Phys. 125, 104109, 2006.
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm.mta;

import au.anu.edu.qm.QMAtom;

import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.ConnectivityBuilder;

public class Fragmentor {
   public def this(rGoodness:Double, maxFragSize:Int, mol:Molecule[QMAtom]) {
       // reorder the atom indices
 
       // next build the connectivity for this molecule
       val conn = new ConnectivityBuilder[QMAtom]();
     
       conn.buildConnectivity(mol); 

       // generate the fragments
   }
}

