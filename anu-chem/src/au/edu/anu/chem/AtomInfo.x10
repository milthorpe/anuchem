/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2012.
 */
package au.edu.anu.chem;

import x10.util.HashMap;

/**
 * This class stores basic information on Atom types.
 * All values are in a.u.
 *
 * @author milthorpe
 */
public class AtomInfo {

    public static class AtomProps(
        symbol:String,
        atomicNumber:Int, 
        covalentRadius:Double, 
        vdWRadius:Double, 
        atomicMass:Double
    ) {

        public def this(symbol:String, atomicNumber:Int, covalentRadius:Double, vdWRadius:Double, atomicMass:Double) {
            property(symbol, atomicNumber, covalentRadius, vdWRadius, atomicMass);
        }
    }

    private static val _instance = new AtomInfo();
    public static def getInstance() {
        return _instance;
    }

    private val bySymbol = new HashMap[String, AtomProps]();

    /** Atom properties stored in order of atomic number.  0 is a placeholder! */
    private val bySpecies:Rail[AtomProps];

    private def this() {
        bySymbol.put("H",  new AtomProps("H",  1,  0.718, 2.268, 1.0079));
        bySymbol.put("He", new AtomProps("He", 2,  0.605, 2.646, 4.0026));
        bySymbol.put("Li", new AtomProps("Li", 3,  2.532, 3.439, 6.941 ));
        bySymbol.put("Be", new AtomProps("Be", 4,  1.701, 2.891, 9.0122));
        bySymbol.put("B",  new AtomProps("B",  5,  1.550, 3.628, 10.81 ));
        bySymbol.put("C",  new AtomProps("C",  6,  1.455, 3.213, 12.011));
        bySymbol.put("N",  new AtomProps("N",  7,  1.417, 2.929, 14.007));
        bySymbol.put("O",  new AtomProps("O",  8,  1.380, 2.872, 15.999));
        bySymbol.put("F",  new AtomProps("F",  9,  1.342, 2.778, 18.998));
        bySymbol.put("Ne", new AtomProps("Ne", 10, 1.304, 2.910, 20.1797));
        bySymbol.put("Na", new AtomProps("Na", 11, 2.910, 4.290, 22.9898));
        bySymbol.put("Mg", new AtomProps("Mg", 12, 2.457, 3.269, 24.3050));
        bySymbol.put("Al", new AtomProps("Al", 13, 2.230, 3.477, 26.9815));
        bySymbol.put("Si", new AtomProps("Si", 14, 2.098, 3.969, 28.085));
        bySymbol.put("P",  new AtomProps("P",  15, 2.003, 3.402, 30.9738));
        bySymbol.put("S",  new AtomProps("S",  16, 2.003, 3.402, 32.06));
        bySymbol.put("Cl", new AtomProps("Cl", 17, 1.871, 3.307, 35.45));
        bySymbol.put("Ar", new AtomProps("Ar", 18, 1.833, 3.553, 39.948));
        bySymbol.put("Br", new AtomProps("Br", 35, 2.154, 3.496, 79.904));

        bySpecies = new Rail[AtomProps](36);
        for (entry in bySymbol.entries()) {
            val props = entry.getValue();
            bySpecies(props.atomicNumber) = props;
        }
    }

    public def getSpecies(speciesId:Int) {
        val species = _instance.bySpecies(speciesId);
        if (species == null) {
            throw new IllegalArgumentException("unknown species " + speciesId);
        }
        return species;
    }

    public def getSpecies(symbol:String) {
        val species = _instance.bySymbol.getOrElse(symbol, null);
        if (species == null) {
            throw new IllegalArgumentException("unknown species " + symbol);
        }
        return species;
    }

    public def getAtomicNumber(symbol:String) = getSpecies(symbol).atomicNumber;
    public def getCovalentRadius(symbol:String) = getSpecies(symbol).covalentRadius;
    public def getVdWRadius(symbol:String) = getSpecies(symbol).vdWRadius;
    public def getAtomicMass(symbol:String) = getSpecies(symbol).atomicMass;

    public def getAtomicNumber(species:Int) = getSpecies(species).atomicNumber;
    public def getCovalentRadius(species:Int) = getSpecies(species).covalentRadius;
    public def getVdWRadius(species:Int) = getSpecies(species).vdWRadius;
    public def getAtomicMass(species:Int) = getSpecies(species).atomicMass;
}

