/* 
 * PhyloTree.java
 *
 * Defines a phylogenetic tree, which is a strictly binary tree 
 * that represents inferred hierarchical relationships between species
 * 
 * There are weights along each edge; the weight from parent to left child
 * is the same as parent to right child.
 *
 * Students may only use functionality provided in the packages
       
 *     
 * Use of any additional Java Class Library components is not permitted 
 * 
 * Audrey Clark and Sameen Ahmad
 *
 */
 
import java.lang.*;
import java.util.*;
import java.io.*;

public class PhyloTree {
   private static PhyloTreeNode overallRoot;    // The actual root of the overall tree
   private int printingDepth;            // How many spaces to indent the deepest 
                                          // node when printing
   private int numSpecies = 0;
   private static ArrayList<Species> descendents = new ArrayList<Species>();

    // CONSTRUCTOR

    // PhyloTree
    // Pre-conditions:
    //        - speciesFile contains the path of a valid FASTA input file
    //        - printingDepth is a positive number
    // Post-conditions:
    //        - this.printingDepth has been set to printingDepth
    //        - A linked tree structure representing the inferred hierarchical
    //          species relationship has been created, and overallRoot points to
    //          the root of this tree
    // Notes:
    //        - A lot happens in this step!  See assignment description for details
    //          on the input format file and how to construct the tree
    //        - If you encounter a FileNotFoundException, print to standard error
    //          "Error: Unable to open file " + speciesFilename
    //          and exit with status (return code) 1
    //        - Most of this should be accomplished by calls to loadSpeciesFile and buildTree
   public PhyloTree(String speciesFile, int printingDepth) {
        //test to see if file is valid, then send to loadSpeciesFile
      Scanner input;
      try{
         //open Scanner
         input = new Scanner(new File(speciesFile));
         //call loadSpeciesFile
         this.numSpecies = 0;
         while(input.hasNext()){
            if(input.next().contains(">")){
               numSpecies ++;
               
            }
         }
         //System.out.println("counting..." + numSpecies);
         Species[] speciesObj = (loadSpeciesFile(speciesFile));
         buildTree(speciesObj);
         descendents = new ArrayList<Species>(Arrays.asList(speciesObj));
      }
      catch(FileNotFoundException ex){
         System.out.println("*Error: Unable to open file " + speciesFile);
         System.exit(1);        
      }
   
      this.printingDepth = printingDepth;
      return;
   }

    // ACCESSORS

    // getOverallRoot
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the overall root
   public PhyloTreeNode getOverallRoot() {
      return this.overallRoot;
   }

    // toString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
    //        (it will be a kind of reverse in-order [RNL] traversal)
    //    - Can be a simple wrapper around the following toString
    //    - Hint: StringBuilder is much faster than repeated concatenation
   public String toString() {
      return toString(this.getOverallRoot(), this.getWeightedHeight(), this.getWeightedHeight());
   }

    // toString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    //    - weightedDepth is the sum of the edge weights from the
    //      overall root to the current root
    //    - maxDepth is the weighted depth of the overall tree
    // Post-conditions:
    //    - Returns a string representation of the tree
    // Notes:
    //    - See assignment description for proper format
   private String toString(PhyloTreeNode node, double weightedDepth, double maxDepth) {
      StringBuilder sb = new StringBuilder();
      String dist = "";
      if(node != null){
         //double k = this.printingDepth*(this.weightedNodeHeight(node)/this.getWeightedHeight());
         int k = (int)(this.printingDepth*(weightedDepth/maxDepth));
         //System.out.println("K: " + k);
         if(!node.isLeaf()){
            //toString(node.rightChild, weightedDepth, maxDepth)
            sb.append(toString(node.getRightChild(), weightedNodeDepth(node.getRightChild()), maxDepth));
         
            //print k periods, k = printingdepth * (this.getWeightedHeight()/this.getHeight())
            for(int i=0; i<=k; i++){
               sb.append(".");
            }
            dist = String.format("%.2f", node.getDistanceToChild());
            sb.append("[NONTERM " + dist + "]" + "\n");
            
            //toString(node.leftChild, weightedDepth, maxDepth)
            sb.append(toString(node.getLeftChild(), weightedNodeDepth(node.getLeftChild()), maxDepth));
         }
         else{
            for(int i=0; i<k; i++){
               sb.append(".");
            }
            sb.append(node.getLabel() + "\n");
         }
      
      }
      return sb.toString();
   }

    // toTreeString 
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for format details
    //    - Can be a simple wrapper around the following toTreeString
   public String toTreeString() {
      return toTreeString(this.getOverallRoot());
   }

    // toTreeString 
    // Pre-conditions:
    //    - node points to the root of a tree you intend to print
    // Post-conditions:
    //    - Returns a string representation in tree format
    // Notes:
    //    - See assignment description for proper format
   private String toTreeString(PhyloTreeNode node) {
      StringBuilder sb = new StringBuilder();
      String val = "";
      if(node.isLeaf()){
         //Print label:weight
         val = String.format("%.5f", node.getParent().getDistanceToChild());
      
         sb.append(node.getLabel() + ":" + val);
         //weight is distance to parent, weight is up to five decimal places
      }
      else{
         sb.append("(");
         sb.append(toTreeString(node.getRightChild()));
         sb.append(",");
         sb.append(toTreeString(node.getLeftChild()));
         if(node == this.getOverallRoot()){
            sb.append(")");
         }
         else{
            val = String.format("%.5f", node.getDistanceToChild());
            sb.append("):" + val);
         }
      }
      return sb.toString();
   }

    // getHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the tree height as defined in class
    // Notes:
    //    - Can be a simple wrapper on nodeHeight
   public int getHeight() {
      return nodeHeight(this.getOverallRoot());
   }

    // getWeightedHeight
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the sum of the edge weights along the
    //      "longest" (highest weight) path from the root
    //      to any leaf node.
    // Notes:
    //   - Can be a simple wrapper for weightedNodeHeight
   public double getWeightedHeight() {
      return weightedNodeHeight(this.getOverallRoot());
   }

    // countAllSpecies
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns the number of species in the tree
    // Notes:
    //    - Non-terminals do not represent species
    //    - This functionality is provided for you elsewhere
    //      just call the appropriate method
   public int countAllSpecies() {
      //Number of leafs = number of non-terminal species?
      return this.getOverallRoot().getNumLeafs();
   }

    // getAllSpecies
    // Pre-conditions:
    //    - None
    // Post-conditions:
    //    - Returns an ArrayList containing all species in the tree
    // Notes:
    //    - Non-terminals do not represent species
   public java.util.ArrayList<Species> getAllSpecies() {
      //getAllDescendantSpecies(this.getOverallRoot(), this.descendents);
      //Collections.sort(descendents, new Comparator<Species>());
      return this.descendents;
   }

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - label is the label of a tree node you intend to find
    //    - Assumes labels are unique in the tree
    // Post-conditions:
    //    - If found: returns the PhyloTreeNode with the specified label
    //    - If not found: returns null
   public PhyloTreeNode findTreeNodeByLabel(String label) {
      return findTreeNodeByLabel(this.getOverallRoot(), label);
   }

    // findLeastCommonAncestor
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns null
    //    - If both nodes can be found: returns the PhyloTreeNode of their
    //      common ancestor with the largest depth
    //      Put another way, the least common ancestor of nodes A and B
    //      is the only node in the tree where A is in the left tree
    //      and B is in the right tree (or vice-versa)
    // Notes:
    //    - Can be a wrapper around the static findLeastCommonAncestor
   public PhyloTreeNode findLeastCommonAncestor(String label1, String label2) {
   
      //search for the node associated with the labels
      //call findLeastCommonAncestor that takes in two nodes --> findLeastCommonAncestor(node1, node2)
      return findLeastCommonAncestor(findTreeNodeByLabel(label1),findTreeNodeByLabel(label2));
   }
    
    // findEvolutionaryDistance
    // Pre-conditions:
    //    - label1 and label2 are the labels of two species in the tree
    // Post-conditions:
    //    - If either node cannot be found: returns POSITIVE_INFINITY
    //    - If both nodes can be found: returns the sum of the weights 
    //      along the paths from their least common ancestor to each of
    //      the two nodes
   public double findEvolutionaryDistance(String label1, String label2) {
      PhyloTreeNode node1 = findTreeNodeByLabel(label1);
      PhyloTreeNode node2 = findTreeNodeByLabel(label2);
   
      if(node1 == null || node2 == null){
         return Double.POSITIVE_INFINITY;
      }
      else{
         double sum = 0.0;
      
               //evolutionary distance can be found by traversing from ancestor to node1 and to node2
      //we will return the sum of weights from the ancestor to node1 and the ancestor to node2
      
         while(nodeDepth(node1) != nodeDepth(node2)){
            //sum the weights along this path, add to sum
            if(nodeDepth(node1) > nodeDepth(node2)){
               node1 = node1.getParent();
               sum += node1.getDistanceToChild();
            }
            else{
               node2 = node2.getParent();
               sum += node2.getDistanceToChild();
            }
         }
         while(node2 != node1){
            //sum the weights along this path, add to sum
            node1 = node1.getParent();
            sum += node1.getDistanceToChild();
            node2 = node2.getParent();
            sum += node2.getDistanceToChild();
         }
         return sum;      
      }
   }

    // MODIFIER

    // buildTree
    // Pre-conditions:
    //    - species contains the set of species for which you want to infer
    //      a phylogenetic tree
    // Post-conditions:
    //    - A linked tree structure representing the inferred hierarchical
    //      species relationship has been created, and overallRoot points to
    //      the root of said tree
    // Notes:
    //    - A lot happens in this step!  See assignment description for details
    //      on how to construct the tree.
    //    - Be sure to use the tie-breaking conventions described in the pdf
    //    - Important hint: although the distances are defined recursively, you
    //      do NOT want to implement them recursively, as that would be very inefficient
   private void buildTree(Species[] species) {
      //create a Hashtable<String, PhyloTreeNode> of nodes, called forestMap
      HashMap<String, PhyloTreeNode> forestMap = new HashMap<String, PhyloTreeNode>();     
      //create distance hashMap using MultiKeyMap class, called distanceMap
      MultiKeyMap distanceMap = new MultiKeyMap();
      
      //filling forestMap
      for(int i=0; i<species.length; i++){
         PhyloTreeNode node = new PhyloTreeNode(null, species[i]);
         forestMap.put(species[i].getName(), node);
      }
      
      //filling distanceMap
      for(int i=0; i<species.length; i++){
         for(int j=0; j<species.length; j++){
            distanceMap.put(species[i].getName(), species[j].getName(), Species.distance(species[i], species[j]));
         }
      }
      //create smallDist
      double smallDist;
      //initialize variables
      double otherDist1 = 0.0;
      double otherDist2 = 0.0;
      String sp1 = "";
      String sp2 = "";
      PhyloTreeNode node1 = null;
      PhyloTreeNode node2 = null;
      PhyloTreeNode parentnode = null;
      String parentlabel = "";
      
      while(forestMap.size() > 1){
         smallDist = Double.MAX_VALUE;
         //find smallest distance between two nodes in the forestMap
      
         for(Map.Entry<String, PhyloTreeNode> entry1 : forestMap.entrySet()){
            for(Map.Entry<String, PhyloTreeNode> entry2 : forestMap.entrySet()){ 
               //if smaller, update smallDist, store the corresponding nodes
               if(distanceMap.getVal(entry1.getKey(), entry2.getKey()) != null){
                  if(distanceMap.getVal(entry1.getKey(), entry2.getKey()) < smallDist && !entry1.getKey().equals(entry2.getKey())){
                  //smallDist is the value/distance associated w/ the two nodes --> uses getVal(String1, String2) from MultiKeyMap                  
                     smallDist = distanceMap.getVal(entry1.getKey(), entry2.getKey());
                  //saving species labels and nodes associated  
                     sp1 = entry1.getKey();
                     sp2 = entry2.getKey();
                     node1 = entry1.getValue();
                     node2 = entry2.getValue();
                  }
               }
            }
         } 
      
         //remove two nodes associated w/ labels from forestMap
         forestMap.remove(sp1);
         forestMap.remove(sp2);
         
         //create parent node (use 2nd constructer from PhyloTreeNode class --> PhyloTreeNode(label, parent, leftChild, rightChild, double distanceToChild) )
         
         //determine left child and right child (left child is alphabetically earlier)
         //create parent label --> leftchildlabel + rightchildlabel
         //insert parent to forestMap
         if(sp1.compareTo(sp2) <= 0){
            //entry1 is alphabetically first
            parentlabel = sp1 + "+" + sp2;
            parentnode = new PhyloTreeNode(parentlabel, null, node1, node2, smallDist/2); 
         }
         else{
            //entry2 is alphabetically first
            parentlabel = sp2 + "+" + sp1;
            parentnode = new PhyloTreeNode(parentlabel, null, node2, node1, smallDist/2);      
         }  
         //setting parent nodes 
         node1.setParent(parentnode);
         node2.setParent(parentnode);
         
         
         //iterate over forestMap and fill distanceMap with distances from parent to other nodes
         for(Map.Entry<String, PhyloTreeNode> entry3 : forestMap.entrySet()){
            //find distance between parentnode and other nodes in forestmap<String, Node>
            if(!entry3.getKey().equals(parentlabel)){
               otherDist1 = (((double)node1.getNumLeafs()) / ((double)(node2.getNumLeafs()) + ((double)node1.getNumLeafs()))) * distanceMap.getVal(sp1, entry3.getKey());
               otherDist2 = (((double)node2.getNumLeafs()) / ((double)(node2.getNumLeafs()) + ((double)node1.getNumLeafs()))) * distanceMap.getVal(sp2, entry3.getKey());
               //store distance found in distanceMap<String, String, double>
               distanceMap.put(parentlabel, entry3.getKey(), otherDist1+otherDist2);
            }
         }
         forestMap.put(parentlabel, parentnode);  
      }
      //set overallRoot, the last and only node in the forest leftover, the last parentnode
      this.overallRoot = parentnode;
      return;
   }
   
   //General descriptions for our reference:    
      //Depth of a node n = length of path from root to n
         //-->number of edges from node to tree's root node
      //Height of a node n = length of longest path from n to a leaf
         //-->number of edges on the longest path between the node and a leaf
      
      //Depth of a tree = depth of deepest node (largest n depth value)
      //Height of a tree = height of root, length from root to any leaf



    // STATIC
    
    // nodeDepth
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the depth of the node within the overall tree
   public static int nodeDepth(PhyloTreeNode node) {
      if(node == null){
         return -1;
      }
      else{
         if(node.getParent() == null){
            return 0;
         }
         else{
            return 1 + nodeDepth(node.getParent());
         }
      }
   }

    // nodeHeight
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns -1
    //    - Else: returns the height subtree rooted at node
   public static int nodeHeight(PhyloTreeNode node) {
      if (node == null) {
         return -1;
      } 
      else {
         if (node.isLeaf()) {
            return 0;
         } 
         else {
            return 1 + Math.max(nodeHeight(node.getLeftChild()), nodeHeight(node.getRightChild()));
         }
      }
   }
   
   //helper function added for benefit of toString()
   //if node == null, return negative infinity
   //if node has parent, continue to find weightednodeDepth
   //else, return 0.0
   public static double weightedNodeDepth(PhyloTreeNode node){
      if(node == null){
         return Double.NEGATIVE_INFINITY;
      }
      if(node.getParent() != null){
         return node.getParent().getDistanceToChild() + weightedNodeDepth(node.getParent());
      }
      else{
         return 0.0;
      }
   }




    // weightedNodeHeight 
    // Pre-conditions:
    //    - node is null or the root of tree (possibly subtree)
    // Post-conditions:
    //    - If null: returns NEGATIVE_INFINITY
    //    - Else: returns the weighted height subtree rooted at node
    //     (i.e. the sum of the largest weight path from node
    //     to a leaf; this might NOT be the same as the sum of the weights
    //     along the longest path from the node to a leaf)
   public static double weightedNodeHeight(PhyloTreeNode node) {
      if(node != null){
         if(!node.isLeaf()){
            double rightH = node.getDistanceToChild() + weightedNodeHeight(node.getRightChild());
            double leftH = node.getDistanceToChild() + weightedNodeHeight(node.getLeftChild());
            if(rightH >= leftH){
               return rightH;
            }
            else{
               return leftH;
            }
         }
         else{
            return 0.0;
         }
      }
      else{
         return Double.NEGATIVE_INFINITY;
      }
   }

    // loadSpeciesFile
    // Pre-conditions:
    //    - filename contains the path of a valid FASTA input file
    // Post-conditions:
    //    - Creates and returns an array of species objects representing
    //      all valid species in the input file
    // Notes:
    //    - Species without names are skipped
    //    - See assignment description for details on the FASTA format
    // Hints:
    //    - Because the bar character ("|") denotes OR, you need to escape it
    //      if you want to use it to split a string, i.e. you can use "\\|" 
   public static Species[] loadSpeciesFile(String filename) {
      String inputline = "";
      int inputlineCount = 0;
      String[] nameArray = new String[7];
      String spname = "";
      String strSeq = "";
      ArrayList<Species> tempArray = new ArrayList<Species>();
      
      try{
         Scanner input = new Scanner(new File(filename));
         inputline = input.next();
         inputlineCount++;
         while(input.hasNext()){          
            nameArray = inputline.split("\\|");
            spname = nameArray[nameArray.length-1];
            inputline = input.next();
            inputlineCount++;
            while(!inputline.contains(">") && input.hasNext()){
                  //Sequence detected
               strSeq = strSeq + inputline;
               inputline = input.next();
               inputlineCount++;
            }
            if((!input.hasNext())){
               strSeq = strSeq + inputline;
            }
               //store sequence
            String[] seqA = new String[strSeq.length()];
            seqA = strSeq.split("");
            strSeq = "";
               //create Species object with name and sequence
            tempArray.add(new Species(spname,seqA));
         }
         
         //System.out.println("inputlineCount: " +inputlineCount);
      }
      catch(FileNotFoundException ex){
         System.out.println("File not found: " + filename);
         System.exit(1);
      }
      
      Species[] speciesObj = tempArray.toArray(new Species[tempArray.size()]);
      descendents = tempArray;      
      String[] strA = new String[speciesObj[5].getSequence().length];
      strA = speciesObj[5].getSequence();
      return speciesObj;
   }
    // getAllDescendantSpecies
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - descendants is a non-null reference variable to an empty arraylist object
    // Post-conditions:
    //    - descendants is populated with all species in the subtree rooted at node
    //      in in-/pre-/post-order (they are equivalent here)
   private static void getAllDescendantSpecies(PhyloTreeNode node, java.util.ArrayList<Species> descendants) {   
      for(int i=0; i<descendents.size(); i++){
         //find node associated with label in Species array
         findTreeNodeByLabel(overallRoot, descendents.get(i).getName());
         //test if node is a leaf
         if(!node.isLeaf()){
            //if node is not a leaf, remove
            descendents.remove(descendents.get(i));
         }
      }
   }

    // findTreeNodeByLabel
    // Pre-conditions:
    //    - node points to a node in a phylogenetic tree structure
    //    - label is the label of a tree node that you intend to locate
    // Post-conditions:
    //    - If no node with the label exists in the subtree, return null
    //    - Else: return the PhyloTreeNode with the specified label 
    // Notes:
    //    - Assumes labels are unique in the tree
   private static PhyloTreeNode findTreeNodeByLabel(PhyloTreeNode node, String label) {
      if(node.getLabel().equals(label)){
         return node;
      }
      if(node.isLeaf()){
         return null;
      }
      else{
         PhyloTreeNode right = findTreeNodeByLabel(node.getRightChild(), label);
         PhyloTreeNode left =findTreeNodeByLabel(node.getLeftChild(), label);
         if(left != null){
            return left;
         }
         else{
            return right;
         }
      }
   }

    // findLeastCommonAncestor
    // Pre-conditions:
    //    - node1 and node2 point to nodes in the phylogenetic tree
    // Post-conditions:
    //    - If node1 or node2 are null, return null
    //    - Else: returns the PhyloTreeNode of their common ancestor 
    //      with the largest depth
   private static PhyloTreeNode findLeastCommonAncestor(PhyloTreeNode node1, PhyloTreeNode node2) {
      if(node1 == null || node2== null){ //if one of the given nodes are null
         return null;
      }
      else{
         //compare depths of node1 and node2
         if(nodeDepth(node1) > nodeDepth(node2)){
            return findLeastCommonAncestor(node1.getParent(), node2);
         }
         else if(nodeDepth(node1) < nodeDepth(node2)){
            return findLeastCommonAncestor(node2.getParent(), node1);
         }
         else if(node1 == node2){
            return node1;
         }
         else{
            return findLeastCommonAncestor(node1.getParent(), node2.getParent());
         
         }
      }
   }
}
