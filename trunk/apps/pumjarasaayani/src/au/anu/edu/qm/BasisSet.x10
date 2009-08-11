/**
 * BasisSet.x10
 *
 * Represents a basis set, used for setting up basis functions (by BasisFunctions)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class BasisSet { 
    var name:String;
    var basisInfo:HashMap[String, AtomicBasis];

    public def this(name:String) { 
       this.name = name;

       init(name);
    } 

    public def this(name:String, basisDir:String) {
       this.name = name;

       init(name, basisDir);
    }

    private def init(name:String, basisDir:String) {
       basisInfo = new HashMap[String, AtomicBasis]();

       // TODO: use file separator!
       val fil = new FileReader(new File(basisDir + "/" + name));       
       fil.close();
    }
 
    private def init(name:String) : void {
       basisInfo = new HashMap[String, AtomicBasis]();

       if (name.equals("sto3g")) {
          val hBasis = new AtomicBasis();
          val hOrb   = new Orbital("S");
          hOrb.add(3.425251, 0.154329);
          hOrb.add(0.623914, 0.535328);
          hOrb.add(0.168855, 0.444635); 
          hBasis.addOrbital(hOrb);
          basisInfo.put("H", hBasis);

          val cBasis = new AtomicBasis();
          val cOrb1  = new Orbital("S");
          cOrb1.add(71.616837, 0.154329);
          cOrb1.add(13.045096, 0.535328);
          cOrb1.add(3.530512, 0.444635);
          cBasis.addOrbital(cOrb1);
          val cOrb2  = new Orbital("S");
          cOrb2.add(2.941249, -0.099967);
          cOrb2.add(0.683483, 0.399513);
          cOrb2.add(0.222290, 0.700115);
          cBasis.addOrbital(cOrb2); 
          val cOrb3  = new Orbital("P");
          cOrb3.add(2.941249, 0.155916);
          cOrb3.add(0.683483, 0.607684);
          cOrb3.add(0.222290, 0.391957);
          cBasis.addOrbital(cOrb3);
          basisInfo.put("C", cBasis);    

          val nBasis = new AtomicBasis();
          val nOrb1  = new Orbital("S");
          nOrb1.add(99.106169, 0.154329);
          nOrb1.add(18.052312, 0.535328);
          nOrb1.add(4.885660, 0.444635);
          nBasis.addOrbital(nOrb1);
          val nOrb2  = new Orbital("S");
          nOrb2.add(3.780456, -0.099967);
          nOrb2.add(0.878497, 0.399513);
          nOrb2.add(0.285714, 0.700115);
          nBasis.addOrbital(nOrb2);
          val nOrb3  = new Orbital("P");
          nOrb3.add(3.780456, 0.155916);
          nOrb3.add(0.878497, 0.607684);
          nOrb3.add(0.285714, 0.391957);
          nBasis.addOrbital(nOrb3);
          basisInfo.put("N", nBasis);    

          val oBasis = new AtomicBasis();
          val oOrb1  = new Orbital("S");
          oOrb1.add(130.709321, 0.154329);
          oOrb1.add(23.808866, 0.535328);
          oOrb1.add(6.443608, 0.444635);
          oBasis.addOrbital(oOrb1);
          val oOrb2  = new Orbital("S");
          oOrb2.add(5.033151, -0.099967);
          oOrb2.add(1.169596, 0.399513);
          oOrb2.add(0.380389, 0.700115);
          oBasis.addOrbital(oOrb2);
          val oOrb3  = new Orbital("P");
          oOrb3.add(5.033151, 0.155916);
          oOrb3.add(1.169596, 0.607684);
          oOrb3.add(0.380389, 0.391957);      
          oBasis.addOrbital(oOrb3);
          basisInfo.put("O", oBasis);    
       } // end if
    }

    public def getName():String = this.name;

    public def getBasis(atom:Atom) : AtomicBasis {
        return (basisInfo.get(atom.getSymbol()) as AtomicBasis);
    }
}

