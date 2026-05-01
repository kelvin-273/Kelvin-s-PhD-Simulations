include "globals.mzn";

int: nLoci;
int: nParents;
set of int: NLoci = 1..nLoci;
set of int: NParents = 1..nParents;
set of int: NGametes = 0..2^nLoci-1;
array[NParents, 1..2] of NGametes: parents;
array[NGametes, NGametes] of var bool: genotypes;

% source nodes must be drawn from parents
% NOTE: this is currently wrong as this is forcing all parents to be used
constraint forall(i in NParents) (genotypes[parents[i, 1], parents[i, 2]]);

% target is created
constraint genotypes[2^nLoci-1, 2^nLoci-1];

% create the potential edges for the dag and dconnected constrants
int: nGenotypes = 2^(2*nLoci);
set of int: NGenotypes = 0..nGenotypes-1;
array[NGenotypes] of var bool: ns = array1d(NGenotypes, genotypes);
array[NGenotypes, NGenotypes] of var bool: left_es;
array[NGenotypes, NGenotypes] of var bool: right_es;
list of NGenotypes: left_from = [ x | x, z in NGenotypes where
  let { NGametes: gx = decodeGenotypeL(z) } in recombinationBoth(x, gx)
];
list of NGenotypes: left_to = [ z | x, z in NGenotypes where
  let { NGametes: gx = decodeGenotypeL(z) } in recombinationBoth(x, gx)
];
list of NGenotypes: right_from = [ x | x, z in NGenotypes where
  let { NGametes: gy = decodeGenotypeR(z) } in recombinationBoth(x, gy)
];
list of NGenotypes: right_to = [ z | x, z in NGenotypes where
  let { NGametes: gy = decodeGenotypeR(z) } in recombinationBoth(x, gy)
];


predicate recombinationBoth(NGenotypes: z, NGametes: gz) = (
  let {
    NGametes: gx = decodeGenotypeL(z),
    NGametes: gy = decodeGenotypeR(z),
  } in recombinationForward(gx, gy, gz) \/ recombinationForward(gy, gx, gz)
);
predicate recombinationForward(NGametes: gx, NGametes: gy, NGametes: gz) = (
  exists (k in NLoci) (
    let {
      int: a1 = 2^k,
      int: a2 = gy mod a1,
      int: a3 = gx div a1,
      int: a4 = a3 * a1,
    } in a2 + a4 = gz
  )
);
function NGenotypes: encodeGenotype(NGametes: gx, NGametes: gy) = (
  2^nLoci * gx + gy
);
function NGametes: decodeGenotypeL(NGenotypes: x) = (
  x div 2^nLoci
);
function NGametes: decodeGenotypeR(NGenotypes: x) = (
  x mod 2^nLoci
);

constraint dconnected(
  left_from, left_to, ns, left_es
);
