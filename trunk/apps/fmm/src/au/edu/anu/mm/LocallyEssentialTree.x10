package au.edu.anu.mm;

/**
 * This class represents the Locally Essential Tree (LET) of
 * a single place.  This is the combined interaction lists
 * of all the boxes assigned to that place.
 * @author milthorpe
 */
public class LocallyEssentialTree {
    //public val combinedUList : ValRail[Point(3)];
    public val combinedVList : ValRail[ValRail[Point(3)]];
    //public val uListMin : ValRail[Int](3);
    //public val uListMax : ValRail[Int](3);
    public val vListMin : ValRail[ValRail[Int](3)];
    public val vListMax : ValRail[ValRail[Int](3)];

    /**
     * A cache of multipole copies for the combined V-list of all
     * boxes at this place.  Used to overlap fetching of the multipole
     * expansions with other computation.
     * The ValRail has one element for each level; each element
     * holds the portion of the combined V-list for that level.
     */
    public global val multipoleCopies : ValRail[Array[Future[MultipoleExpansion]](3)];
    
    public def this(//combinedUList : ValRail[Point(3)],
                combinedVList : ValRail[ValRail[Point(3)]],
                //uListMin : ValRail[Int](3),
                //uListMax : ValRail[Int](3),
                vListMin : ValRail[ValRail[Int](3)],
                vListMax : ValRail[ValRail[Int](3)]) {
        //this.combinedUList = combinedUList;
        this.combinedVList = combinedVList;
        //this.uListMin = uListMin;
        //this.uListMax = uListMax;
        this.vListMin = vListMin;
        this.vListMax = vListMax;
        val multipoleCopies = Rail.make[Array[Future[MultipoleExpansion]](3)](combinedVList.length());
        for ((i) in 0..combinedVList.length()-1) {
            val multipoleCopiesLevelRegion : Region(3) = [vListMin(i)(0)..vListMax(i)(0), vListMin(i)(1)..vListMax(i)(1), vListMin(i)(2)..vListMax(i)(2)];
            multipoleCopies(i) = Array.make[Future[MultipoleExpansion]](multipoleCopiesLevelRegion);
        }
        this.multipoleCopies = multipoleCopies;
    }
    
}
