/*
3D kd tree for nearest neighbour finding
 Based on examples from https://rosettacode.org/wiki/K-d_tree#D

@author: Nick Gibbons
*/
module geom.misc.kdtree;

import std.stdio, std.algorithm, std.math, std.random;

struct Node {
    double[3] x;
    Node* left, right;
    size_t blkid, cellid;
}

void quicksort(size_t idx)(Node[] nodes) pure nothrow @nogc {
/*
   Sketchy quicksort based on C implementation from:
      https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#D
   The goal is to sort an array of nodes by one of their dimensions.
   For example, quicksort!(0)(nodes) sorts by the x direction, index 0
   quicksort!(1)(nodes) sorts by the y direction, index 1

   @author: Nick Gibbons
*/
  size_t len = nodes.length;
  if (len < 2) return;

  immutable double pivot = nodes[len / 2].x[idx];

  size_t i, j;
  for (i = 0, j = len - 1; ; i++, j--) {
      while (nodes[i].x[idx] < pivot) i++;
      while (nodes[j].x[idx] > pivot) j--;

      if (i >= j) break;

      swap(nodes[i].x, nodes[j].x);
      swap(nodes[i].blkid, nodes[j].blkid);
      swap(nodes[i].cellid, nodes[j].cellid);
  }

  quicksort!(idx)(nodes[0 .. i]);
  quicksort!(idx)(nodes[i .. $]);
}

Node* makeTree(size_t i = 0)(Node[] nodes) pure nothrow @nogc {
    if (nodes.length==0) return null;

    quicksort!(i)(nodes);
    size_t medianIdx = nodes.length/2;
    Node* n = &nodes[medianIdx];

    // medianIdx will be zero if there is only element in nodes
    // In that case we can just skip this bit and return n
    if (medianIdx>0) {
        enum i2 = (i + 1) % 3;
        immutable size_t nPos = n - nodes.ptr;
        n.left = makeTree!(i2)(nodes[0 .. medianIdx]);
        n.right = makeTree!(i2)(nodes[medianIdx + 1 .. $]);
    }

    return n;
}

double distance_squared(in Node a, in Node b) pure nothrow @nogc {
        double result = (a.x[0] - b.x[0]) ^^ 2
                      + (a.x[1] - b.x[1]) ^^ 2
                      + (a.x[2] - b.x[2]) ^^ 2;
        return result;
    }

void fast_nearest(in Node* root,
                  in Node nd,
                  in size_t i,
                  ref const(Node)* best,
                  ref double bestDist,
                  ref size_t nVisited) pure nothrow @nogc {
/*
   Fast searching for nearest neighbour in 3D, modified from
      https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#D
   Note that this routine returns "bestDist" the closest distance, SQUARED
*/

    if (root == null)
        return;

    immutable double d = distance_squared(*root, nd);
    immutable double dx = root.x[i] - nd.x[i];
    immutable double dx2 = dx ^^ 2;
    nVisited++;

    if (!best || d < bestDist) {
        bestDist = d;
        best = root;
    }

    // If chance of exact match is high.
    if (!bestDist)
        return;

    immutable i2 = (i + 1 >= 3) ? 0 : i + 1;

    fast_nearest(dx > 0 ? root.left : root.right,
                 nd, i2, best, bestDist, nVisited);
    if (dx2 >= bestDist)
        return;
    fast_nearest(dx > 0 ? root.right : root.left,
                nd, i2, best, bestDist, nVisited);
}

double slow_nearest(Node[] nodes, Node point){
/*
   Check the distance to the nearest point by literally checking all of them
*/

    double cx,cy,cz,dsq;
    double dsqmin = 1e32;
    foreach(node; nodes){
        cx = (node.x[0] - point.x[0]);
        cy = (node.x[1] - point.x[1]);
        cz = (node.x[2] - point.x[2]);
        dsq = cx*cx + cy*cy + cz*cz;
        dsqmin = fmin(dsq, dsqmin);
    }
    return sqrt(dsqmin);
}

/// test sorting
unittest {
    Node[] nodes = [{[2, 3, 0]}, {[5, 4, 0]}, {[9, 6, 0]},
                    {[4, 7, 0]}, {[8, 1, 0]}, {[7, 2, 0]}];

    Node[] targetnodesx= [{[2, 3, 0]}, {[4, 7, 0]}, {[5, 4, 0]},
                         {[7, 2, 0]}, {[8, 1, 0]}, {[9, 6, 0]}];
    Node[] targetnodesy= [{[8, 1, 0]}, {[7, 2, 0]}, {[2, 3, 0]},
                          {[5, 4, 0]}, {[9, 6, 0]}, {[4, 7, 0]}];
    quicksort!(0)(nodes);
    assert(nodes==targetnodesx);
    quicksort!(1)(nodes);
    assert(nodes==targetnodesy);
}

/// test kdtree construction
unittest {
    Node[] nodes = [{[2, 3, 1]}, {[5, 4, 2]}, {[9, 6, 3]},
                    {[4, 7, 6]}, {[8, 1, 8]}, {[7, 2, 4]}];

    auto root = makeTree(nodes);
    //writeln("            ", root.x);
    //writeln("     ",root.left.x, root.right.x);
    //writeln(root.left.left.x, root.left.right.x, root.right.left.x, root.right.right);
    assert(root.left.x[0] <= root.x[0]);
    assert(root.x[0] <= root.right.x[0]);
    assert(root.left.left.x[1] <= root.left.x[1]);
    assert(root.left.right.x[1] >= root.left.x[1]);
    assert(root.right.left.x[1] <= root.right.x[1]);
}

// test nearest computation
unittest {
    Node[] nodes = [{[0.0, 4.0, 0.0]}, {[1.0, 4.0, 0.0]}, {[2.0, 4.0, 0.0]}, {[3.0, 4.0, 0.0]}, {[4.0, 4.0, 0.0]},
                     {[0.0, 4.0, 2.0]}, {[1.0, 4.0, 2.0]}, {[2.0, 4.0, 2.0]}, {[3.0, 4.0, 2.0]}, {[4.0, 4.0, 2.0]},
                     {[0.0, 4.0, 3.0]}, {[1.0, 4.0, 3.0]}, {[2.0, 4.0, 3.0]}, {[3.0, 4.0, 3.0]}, {[4.0, 4.0, 3.0]},
                     {[0.0, 4.0, 4.0]}, {[1.0, 4.0, 4.0]}, {[2.0, 4.0, 4.0]}, {[3.0, 4.0, 4.0]}, {[4.0, 4.0, 4.0]},
                     {[0.0, 4.0, 1.0]}, {[1.0, 4.0, 1.0]}, {[2.0, 4.0, 1.0]}, {[3.0, 4.0, 1.0]}, {[4.0, 4.0, 1.0]}];

    foreach(i, ref node; nodes) {
        node.blkid = 0;
        node.cellid = i;
    }

    auto root = makeTree(nodes);
    size_t ntests=100;
    auto rng = Random(19920829);
    double rngx,rngy,rngz;

    foreach(n; 0 .. ntests){
        rngx = uniform(0.0, 4.0, rng);
        rngy = uniform(0.0, 4.0, rng);
        rngz = uniform(0.0, 4.0, rng);
        Node thisPt = {[rngx, rngy, rngz]};
        double actual_distance = slow_nearest(nodes, thisPt);

        const(Node)* found = null;
        double bestDist = 0;
        size_t nVisited = 0;
        root.fast_nearest(thisPt, 0, found, bestDist, nVisited);

        writefln("point: %s, nearest: %s, dist = %g actual = %g id = %d",
                 thisPt.x, found.x, sqrt(bestDist), actual_distance, found.cellid);
        assert(isClose(sqrt(bestDist), actual_distance, 1e-12));
    }
}

//int main() {return 0;}
