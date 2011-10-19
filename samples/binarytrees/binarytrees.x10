/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */

/**
 * Binary trees benchmark from:
 * The Computer Language Benchmarks Game
 * http://shootout.alioth.debian.org/
 *
 * This code is adapted from the Java 6 averaged version, 
 * which was contributed by Jarkko Miettinen using
 * a steady state approximation by Isaac Gouy
 *
 */
public class binarytrees {
    public static def main(args: Rail[String]) {
        var n:Int = 0;
        if (args.size > 0) n = Int.parse(args(0));
        for (i in 0..64) binarytrees.run(n,false);
        binarytrees.run(n,true);
    }

    private static minDepth = 4;

    public static def run(n:Int, isWarm:boolean) {
        val maxDepth = (minDepth + 2 > n) ? minDepth + 2 : n;
        val stretchDepth = maxDepth + 1;

        var check:Int = (TreeNode.bottomUpTree(0,stretchDepth)).itemCheck();
        
        if (isWarm)
         Console.OUT.println("stretch tree of depth " + stretchDepth + "\t check: " + check);
      
        val longLivedTree = TreeNode.bottomUpTree(0,maxDepth);
      
        for (var depth:Int=minDepth; depth<=maxDepth; depth+=2){
            val iterations = 1 << (maxDepth - depth + minDepth);
            check = 0;

            for (i in 1..iterations){
                check += (TreeNode.bottomUpTree(i,depth)).itemCheck();
                check += (TreeNode.bottomUpTree(-i,depth)).itemCheck();
            }
            if (isWarm)
                Console.OUT.println((iterations*2) + "\t trees of depth " + depth + "\t check: " + check);
        }
        if (isWarm)
            Console.OUT.println("long lived tree of depth " + maxDepth + "\t check: "+ longLivedTree.itemCheck());
    }
   
   
    private static class TreeNode {
        private val left:TreeNode;
        private val right:TreeNode;
        private val item:Int;

        def this(left:TreeNode, right:TreeNode, item:Int){ 
             this.left = left;
             this.right = right;
             this.item = item;
        }

        def this(item:Int) {
            this(null, null, item);
        }

        private static def bottomUpTree(item:Int, depth:Int):TreeNode {
            if (depth > 0){
                return new TreeNode(
                        bottomUpTree(2*item-1, depth-1),
                        bottomUpTree(2*item, depth-1),
                        item
                );
            } else {
                return new TreeNode(item);
            }
        }

        private def itemCheck():Int {
             if (left==null) return item;
             else return item + left.itemCheck() - right.itemCheck();
        }
    }
}
