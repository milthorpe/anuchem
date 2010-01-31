/**
 * BasisSet.x10
 *
 * Represents a basis set, used for setting up basis functions (by BasisFunctions)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.io.*;
import x10.util.HashMap;
import au.edu.anu.chem.Atom;

public class BasisSet { 
    global var name:String;
    global var basisInfo:HashMap[String, AtomicBasis{self.at(this)}];

    public def this() { }

    public def make(name:String) { 
       this.name = name;

       basisInfo = new HashMap[String, AtomicBasis{self.at(this)}]();
       init(name);
    } 

    public def make(name:String, basisDir:String) {
       x10.io.Console.OUT.println("\tReading in basis info. for " + name + " from " + basisDir);

       basisInfo = new HashMap[String, AtomicBasis{self.at(this)}]();
       this.name = name;

       try {
         init(name, basisDir);
       } catch(e:Exception) {
         x10.io.Console.OUT.println("Unable to read basis from : " + basisDir + ".");
         x10.io.Console.OUT.println("Will use sto3g basis.");

         init("sto3g");
       } // end of try .. catch block
    }

    private def init(name:String, basisDir:String) throws IOException {
       // TODO: use file separator in the following line!
       val fil = new FileReader(new File(basisDir + "/" + name));       
       val noOfAtoms = Int.parseInt(fil.readLine());

       for(var i:Int=0; i<noOfAtoms; i++) {
          // var words:ArrayList[String] = Utility.split(fil.readLine(), ' ');
          val words = fil.readLine().split(" ");
          var symbol:String = words(0);
          var noOfContractions:Int = Int.parseInt(words(1));

          var atomBasis:AtomicBasis{self.at(this)} = new AtomicBasis();
          atomBasis.make();

          for(var j:Int=0; j<noOfContractions; j++) {
             val words1 = fil.readLine().split(" ");
             var orbitalType:String = words1(0);
             var noOfPrimitives:Int = Int.parseInt(words1(1));

             val orbital:Orbital{self.at(this)} = new Orbital();
             orbital.make(orbitalType);

             for(var k:Int=0; k<noOfPrimitives; k++) { 
                val words2 = fil.readLine().split(" ");
                
                orbital.add(Double.parseDouble(words2(0)),
                            Double.parseDouble(words2(1)));
             } // end for

             atomBasis.addOrbital(orbital); 
          } // end for

          basisInfo.put(symbol, atomBasis);
       } // end for

       fil.close();
    }
 
    private def init(name:String) {
       if (name.equals("sto3g")) {
          var orbType:String;
          val hBasis = new AtomicBasis();
          hBasis.make();
          val hOrb   = new Orbital();
          orbType = "S";
          hOrb.make(orbType as String{self.at(this)});
          hOrb.add(3.425251, 0.154329);
          hOrb.add(0.623914, 0.535328);
          hOrb.add(0.168855, 0.444635); 
          hBasis.addOrbital(hOrb);
          basisInfo.put("H", hBasis);

          val cBasis = new AtomicBasis();
          cBasis.make();
          val cOrb1  = new Orbital();
          orbType = "S";
          cOrb1.make(orbType as String{self.at(this)});
          cOrb1.add(71.616837, 0.154329);
          cOrb1.add(13.045096, 0.535328);
          cOrb1.add(3.530512, 0.444635);
          cBasis.addOrbital(cOrb1);
          val cOrb2  = new Orbital();
          orbType = "S";
          cOrb2.make(orbType as String{self.at(this)});
          cOrb2.add(2.941249, -0.099967);
          cOrb2.add(0.683483, 0.399513);
          cOrb2.add(0.222290, 0.700115);
          cBasis.addOrbital(cOrb2); 
          val cOrb3  = new Orbital();
          orbType = "P";
          cOrb3.make(orbType as String{self.at(this)});
          cOrb3.add(2.941249, 0.155916);
          cOrb3.add(0.683483, 0.607684);
          cOrb3.add(0.222290, 0.391957);
          cBasis.addOrbital(cOrb3);
          basisInfo.put("C", cBasis);    

          val nBasis = new AtomicBasis();
          nBasis.make();
          val nOrb1  = new Orbital();
          orbType = "S";
          nOrb1.make(orbType as String{self.at(this)});
          nOrb1.add(99.106169, 0.154329);
          nOrb1.add(18.052312, 0.535328);
          nOrb1.add(4.885660, 0.444635);
          nBasis.addOrbital(nOrb1);
          val nOrb2  = new Orbital();
          orbType = "S";
          nOrb2.make(orbType as String{self.at(this)});
          nOrb2.add(3.780456, -0.099967);
          nOrb2.add(0.878497, 0.399513);
          nOrb2.add(0.285714, 0.700115);
          nBasis.addOrbital(nOrb2);
          val nOrb3  = new Orbital();
          orbType = "P";
          nOrb3.make(orbType as String{self.at(this)});
          nOrb3.add(3.780456, 0.155916);
          nOrb3.add(0.878497, 0.607684);
          nOrb3.add(0.285714, 0.391957);
          nBasis.addOrbital(nOrb3);
          basisInfo.put("N", nBasis);    

          val oBasis = new AtomicBasis();
          oBasis.make();
          val oOrb1  = new Orbital();
          orbType = "S";
          oOrb1.make(orbType as String{self.at(this)});
          oOrb1.add(130.709321, 0.154329);
          oOrb1.add(23.808866, 0.535328);
          oOrb1.add(6.443608, 0.444635);
          oBasis.addOrbital(oOrb1);
          val oOrb2  = new Orbital();
          orbType = "S";
          oOrb2.make(orbType as String{self.at(this)});
          oOrb2.add(5.033151, -0.099967);
          oOrb2.add(1.169596, 0.399513);
          oOrb2.add(0.380389, 0.700115);
          oBasis.addOrbital(oOrb2);
          val oOrb3  = new Orbital();
          orbType = "P";
          oOrb3.make(orbType as String{self.at(this)});
          oOrb3.add(5.033151, 0.155916);
          oOrb3.add(1.169596, 0.607684);
          oOrb3.add(0.380389, 0.391957);      
          oBasis.addOrbital(oOrb3);
          basisInfo.put("O", oBasis);    
       } // end if
    }

    public def getName() = this.name;

    public def getBasis(atom:Atom{self.at(this)}) = basisInfo.getOrElse(atom.symbol, null);
}

