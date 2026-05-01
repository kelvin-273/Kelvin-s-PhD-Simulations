set of int: B = {0,1};

int: nLoci;
set of int: NLoci = 1..nLoci;
set of int: NLociDiff = 1..nLoci-1;

int: nPop;
set of int: NPop = 1..nPop;

array[NPop, NLoci] of B: A;
array[NLociDiff] of float: r;

array[NPop, NLoci] of var B: x;
array[NPop, NLociDiff] of var B: y;
array[NPop] of var B: z;

constraint forall (i in NLoci, j in NPop) (
  x[j, i] <= A[j, i]
);
constraint forall (i in NLoci) (
  sum(x[NPop, i]) = 1
);
constraint forall (i in NLoci, j in NPop) (
  y[j,i] >= x[j,i] - x[j,i+1]
);
constraint forall (j in NPop) (
  M * z[j] >= sum(x[NLoci, j])
);

float: a, b;
var float: rho;
constraint sum(i in NLoci, j in NPop) (

);
