/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.util;

import x10.util.Box;

/**
 * This class implements a scalable tree barrier across all X10 places.
 * TODO: this functionality should be provided by clocks,
 * but they are currently inefficient as they are implemented using
 * point-to-point rather than collective communications.
 * TODO: deadlocks due to use of X10 conditional atomic blocks: see XTENLANG-1660
 * @see  Mellor-Crummey and Scott (1991). Algorithms for Scalable 
 * Synchronization on Shared-Memory Multiprocessors. ACM TOCS, February 1991.
 * http://www.cs.rochester.edu/u/scott/papers/1991_TOCS_synch.pdf 
 */
final public class ScalableTreeBarrier {
    final static class TreeNode {
        var sense : Boolean = true;
        var parentSense : Boolean = false;
        var parent : Box[GlobalRef[TreeNode]];
        val children = new Array[Box[GlobalRef[TreeNode]]](0..1);
        val haveChild = new Array[Boolean](0..3);
        val childNotReady = new Array[Boolean](0..3);

        public def this(numPlaces : Int) {
            this.parent = parent;
            val i = here.id;
            for (j in 0..3) {
                haveChild(j) = (4*i+j+1 < numPlaces);
                childNotReady(j) = haveChild(j);
            }
        }
    }

    val nodes : DistArray[TreeNode](1);

    public def this() {
        this(PlaceGroup.WORLD);
    }

    public def this(places : PlaceGroup) {
        val nodes = DistArray.make[TreeNode](Dist.makeUnique(places), ([i] : Point) => new TreeNode(places.size));
        val P = nodes.region.size();
        finish ateach([i] in nodes) {
            val thisNode = nodes(i);
            val parentId = Math.floor((i-1)/4) as Int;
            thisNode.parent = (i == 0) ? null : at(nodes.dist(parentId)) {new Box[GlobalRef[TreeNode]](GlobalRef[TreeNode](nodes(parentId)))};
            val childId0 = 2*i+1;
            thisNode.children(0) = (childId0 >= P) ? null : at(nodes.dist(childId0)) {new Box[GlobalRef[TreeNode]](GlobalRef[TreeNode](nodes(childId0)))};
            val childId1 = 2*i+2;
            thisNode.children(1) = (childId1 >= P) ? null : at(nodes.dist(childId1)) {new Box[GlobalRef[TreeNode]](GlobalRef[TreeNode](nodes(childId1)))};
        }
        this.nodes = nodes;
    }

    public def barrier() {
        //Console.OUT.println("treeBarrier starting at " + here.id);

        val thisNode = nodes(here.id);
        //Console.OUT.println("waiting for children at " + here.id);

        val booleanOr = (a:Boolean,b:Boolean) => a || b;
        when (thisNode.childNotReady.reduce(booleanOr, true));

        //Console.OUT.println("children ready at " + here.id);
        // prepare for next barrier
        Array.copy(thisNode.haveChild, thisNode.childNotReady);

        val i = here.id;
        if (i != 0) {
            // let parent know I'm ready
            Console.OUT.println("I am ready : " + i);
            val parent = thisNode.parent();
            at(parent) {
                Console.OUT.println("At : " + here.id + " setting " + ((i-1)%4));
                parent().childNotReady((i-1)%4) = false;
                Console.OUT.println("At : " + here.id + " have set " + ((i-1)%4));
            }
            Console.OUT.println("Now sleep : " + i);

            // wait until my parent signals wakeup
            when (thisNode.parentSense == thisNode.sense);

            Console.OUT.println("Woke up : " + i);
        }
        // signal children in wakeup tree
        val thisNodeSense = thisNode.sense;
        for (p in thisNode.children) {
            val childBox = thisNode.children(p);
            if (childBox != null) {
                val child = childBox();
                Console.OUT.println("About to signal : " + child.home);
                at(child) { 
                    child().parentSense = thisNodeSense;
                    Console.OUT.println("Signalled : " + here.id);
                }
            }
        }
        thisNode.sense = !thisNode.sense;
        Console.OUT.println("treeBarrier completed at " + here.id);
    }
}


