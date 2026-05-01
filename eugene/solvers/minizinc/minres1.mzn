include "modelGenotypes.mzn";

set of float: P = 0.0..1.0;
set of int: NLoci0 = NLoci union {0};
array[NLociDiff] of P: recombinationRate;
P: gamma;
constraint forall (j in NLociDiff) (
  0 <= recombinationRate[j] /\ recombinationRate[j] <= 0.5
);

array[NLoci] of P: cumulativeRecRate = [product(j in NLociDiff)(1-recombinationRate[j])] ++ [
  recombinationRate[k] * product(j in NLociDiff where j != k) (1-recombinationRate[j])
  | k in NLociDiff];

/*
array[NTreeCells, 1..2, NLoci] of var 1..2: source;
array[NTreeCells, 1..2, NLociDiff] of var B: swap;
constraint forall (i in NTreeCells, c in 1..2 where treeType[i] != Node) (
  forall (j in NLoci) (source[i, c, j] = 1) /\
  forall (j in NLociDiff) (swap[i, c, j] = 1)
);
constraint forall (i in NTreeCells where treeType[i] == Node) (
  forall (j in NLociDiff) (
    treeData[i, 1, j] == treeData[treeLeft[i], source[i, 1, j], j] /\
    treeData[i, 2, j] == treeData[treeRight[i], source[i, 2, j], j] /\
    treeData[i, 1, j-1] == treeData[treeLeft[i], source[i, 1, j-1], j-1] /\
    treeData[i, 2, j-1] == treeData[treeRight[i], source[i, 2, j-1], j-1] /\
    swap[i, 1, j] = (source[i, 1, j-1] != source[i, 1, j]) /\
    swap[i, 2, j] = (source[i, 2, j-1] != source[i, 2, j]) /\
    true
  )
);
constraint forall (i in NTreeCells where treeType[i] == Node) (
  count(swap[i, 1, NLociDiff], 1) <= 1 /\
  sum(swap[i, 1, NLociDiff]) <= 1
);
*/

function var float: resourcesCrossing(
  var NTreeCells union {0}: igz,
  var NTreeCells union {0}: igx,
  var NTreeCells union {0}: igy,
) = (
  1.0
);
%function var P: probGamete(
%  var NTreeCells: igx,
%  var NTreeCells: igz,
%  1..2: cgz,
%) = (
%  % TODO: replace this way of computing the probabilities with the SPC method <17-01-25>
%  product(j in NLociDiff) (
%    recombinationRate[j] * swap[igz, cgz, j] +
%    (1 - recombinationRate[j]) * (1 - swap[igz, cgz, j])
%  )
%);

% Lengths of a prefix/suffix of x_1/2 that matches that of g^x.
% E.g. k1pref is the length of a prefix of x_1 that matches g^x.
% In maximising the objective, these values will be the lengths of the largest matching prefix/suffix.
array[NTreeCells, 1..2] of var NLoci0: lenPref;
array[NTreeCells, 1..2] of var NLoci0: lenSuff;
constraint forall (i in NTreeCells where treeType[i] == Node) (
  forall (c in 1..2, j in NLoci) (
    (treeData[i, c, j] != treeData[i, treeLeft[i], j] -> lenPref[i, c] < j) /\
    (treeData[i, c, j] != treeData[i, treeLeft[i], j] -> lenSuff[i, c] <= nLoci - j) /\
    true
  )
);
constraint symmetry_breaking_constraint(
  forall (i in NTreeCells, c in 1..2) (
    lenPref[i, c] + lenSuff[i, c] >= nLoci ->
    lenPref[i, c] = nLoci /\ lenSuff[i, c] = nLoci
));

array[NLoci] of P: probRecomb::output = [
  product (j in NLociDiff) (
    if j == k
    then recombinationRate[j]
    else (1-recombinationRate[j])
    endif
  ) | k in NLociDiff] ++ 
  [product(j in NLociDiff) (1-recombinationRate[j])];

/*
Given a starting chromosome c, probRange[lenPref, lenSuff] gives the probability of (x -> g^x)
*/
array[NLoci0, NLoci0] of 0.0..1.0: probRange::output = array2d(NLoci0, NLoci0,
  [
    if a == nLoci /\ b == nLoci
    then 1.0
    elseif a + b < nLoci 
    then 0.0
    elseif a == nLoci
    then sum(k in nLoci-b..nLoci) (probRecomb[k])
    elseif b == nLoci
    then sum(k in 1..a) (probRecomb[k])
    else sum(k in nLoci-b..a) (probRecomb[k]) endif
  | a in NLoci0, b in NLoci0]
);

set of int: NLoci0Squared = 0..(nLoci+1)^2-1;
function 0..(nLoci+1)^2-1: key(NLoci0: a, NLoci0: b) = ((nLoci + 1)*a + b);
any: resourcesLookup = [
  (
    key(x1pref, x1suff), key(x2pref, x2suff),
    key(y1pref, y1suff), key(y2pref, y2suff)
  ): let {
    0.0..1.0: px = 1.0;
    0.0..1.0: py = 1.0;
  } in if px == 1.0 /\ py == 1.0 then 1.0
    elseif gamma == 1.0 \/ px == 0.0 \/ py == 0.0 then infinity
    elseif gamma == 0.0 then 1.0
    else ceil(log2(1-gamma) / log2(1-px*py))
    endif
  | x1pref, x1suff, x2pref, x2suff, y1pref, y1suff, y2pref, y2suff in NLoci0
];


function var 0.0..1.0: probGamete(
  var NTreeCells: igx,
  var NTreeCells: igz,
  1..2: cgz,
) = (
  1.0
  %(probRange[lenPref[igz, 1], lenSuff[igz, 2]] + probRange[lenPref[igz, 2], lenSuff[igz, 1]]) / 2
);

/*
 * Objective
 */

var 0.0..1000.0: obj::output = sum(i in NTreeCells) (
  if treeType[i] == Node
  then resourcesCrossing(i, treeLeft[i], treeRight[i])
  else 0.0
  endif
);
%solve satisfy;
solve minimize obj;
