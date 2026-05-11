(* ::Package:: *)

(* ::Title:: *)
(*Test  Data  for  FusionRings  Package*)


<<Anyonica`
codes = Import["/home/gert/Projects/FusionRings/test/testdata/test_fusion_ring_codes.wdx"];


(* Exporting data *)
tojson[ input_, output_, info_ ] :=
	<|
		"Input" -> input, 
		"Output" -> output,
		"Info" -> info
	|>;

exportjson[ fn_, dt_ ] := 
	Export["/home/gert/Projects/FusionRings/test/testdata/"<>fn<>".json",dt,"JSON"];


(* ::Section::Closed:: *)
(*Creation*)


(* ::Text:: *)
(*Check that creation functions of FusionRings.jl give same result as those of Anyonica*)


zntables = 
	tojson[
		ToString/@Range[8],
		MT[FusionRingZn[#]]&/@Range[8],
		"Input: int i, Output: multiplication table for ring of cyclotomic group of order i"
	];
exportjson["zn_tables",zntables]


smallgroups = 
	Join[Flatten[
		Table[ {"AbelianGroup",{i,j}},{i,2,3},{j,2,3} ],1],
		Table[{"DihedralGroup",i}, {i,3,5 }],
		Table[{"SymmetricGroup",i}, {i,4 }]
	];
smallgroupsmulttabs = FiniteGroupData[#,"MultiplicationTable"]&/@smallgroups;

grouptables = 
	tojson[
		smallgroupsmulttabs,
		MT@*FusionRingFromGroup/@smallgroupsmulttabs,
		"Input: mult tab of small group, Output: mult tab of fusion ring of that group"
	];
exportjson["grouptables",grouptables];


smallabeliangroups = Flatten[Table[ FiniteGroupData[{"AbelianGroup",{i,j}},"MultiplicationTable"],{i,2,3},{j,2,3} ],1];
hitables = 
	tojson[
		smallabeliangroups,
		MT@*FusionRingHI/@smallabeliangroups,
		"Input: mult tab of small abelian group, Output: mult tab of HI fusion ring of that group"
	];
exportjson["HI_tables",hitables]


tytables =
	tojson[
		smallgroupsmulttabs,
		MT@*FusionRingTY/@smallgroupsmulttabs,
		"Input: mult tab of small group, Output: mult tab of TY fusion ring of that group"
	];
exportjson["TY_tables",tytables];


psu2ktables = 
	tojson[
		Range[8],
		MT@*FusionRingPSU2k /@ Range[8],
		"Input: int i, Output: mult tab of PSU(2)_i fusion ring"
	];
	
exportjson["psu2k_tables",psu2ktables];


su2ktables = 
	tojson[
		Range[8],
		MT@*FusionRingSU2k /@ Range[8],
		"Input: int i, Output: mult tab of SU(2)_i fusion ring"
	];
	
exportjson["su2k_tables",su2ktables];


son2tables = 
	tojson[
		Range[8],
		MT@*FusionRingSON2 /@ Range[8],
		"Input: int i, Output: mult tab of SO(i)_2 fusion ring"
	];
exportjson["son2_tables",son2tables];





(* ::Section:: *)
(*Operations*)


(* ::Subsection:: *)
(*Permutations*)


(* Generate 2 random permutations for 64 rings with rank > 3 *)
minrank4codes = Cases[ codes, c_ /; c[[1]] > 3 ];

{ multfreecodes, multcodes } = BinSplit[ minrank4codes, #[[2]] == 1& ];

SeedRandom[1];

randomcodes = Join[ RandomSample[ multfreecodes, 32 ], RandomSample[ multcodes, 32 ] ];

randomperm[ code_ ] := 
	Prepend[1] /@ 
	RandomSample[ 
		Permutations[ Rest @ Range @ First @ code ],
		2
	];
	
perms = randomperm /@ randomcodes;

input = Transpose @ { randomcodes, perms }; 

output = 
	Table[ 
		MT[PermutedRing[FusionRing[ "MultiplicationTable"->inp[[1]]],#]]&/@inp[[2]],
		{inp,MapAt[ MT@*FRBC, input, {All,1} ]}
	];

permtabs = 
	tojson[ 
		input,
		output,
		"Input: list { code, perms } where code is formal code of fusion ring stored in package and perms is list of 2 permutations, Output: list of mult tabs corresponding to permutations of the stored ring"
	];
	
exportjson[ "permuted_tabs", permtabs ]


(* ::Subsection:: *)
(*Tensor products*)


(* tensor products between rings. We can only test for small rings *)
maxrank4codes = Cases[ codes, c_ /; c[[1]] < 5 ];

{ multfreecodes, multcodes } = BinSplit[ maxrank4codes, #[[2]] == 1& ];

SeedRandom[1];

randomcodes = Join[ RandomSample[ multfreecodes, 4 ], RandomSample[ multcodes, 4 ] ];

input = 
	Map[ FRBC, Tuples[{ randomcodes, randomcodes }], {2} ];

output = PMap[ MT[TensorProduct@@#]&, input ];

tptabs = tojson[ 
	Map[ FC, input, {2} ],
	output,
	"Input: couple of formal codes of 2 stored fusion rings, Output: mult tab of tensor product of the rings"
];

exportjson[ "tensor_product_tables", tptabs ]


(* ::Section:: *)
(*Properties*)


(* ::Subsection:: *)
(*sub rings*)


input = codes;
output = SubFusionRings@*FRBC/@codes/.r_FusionRing:>FC[r];

subrings = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, Output: either (1) empty list if ring has no non-trivial sub fusion rings or (2) list of couples ( els, fc ) where els are the elements in the parent ring that make up the subring and fc is the formal code of the subring"
];

exportjson[ "sub_ring_data", subrings ]


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = Sort@*Rest@*FRA/@inputrings;

autos = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, Output: non-trivial fusion ring automorphisms as permutation vectors, sorted lexicographically"	
];

exportjson[ "automorphisms", autos ]


(* ::Subsection:: *)
(*number non zero struct constants*)


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = NNZSC/@inputrings;

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, Output: number non-zero structure constants"
];

exportjson[ "num_nonzero_structure_constants", json ]


(* ::Subsection:: *)
(*FPdims (numerical test)*)


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = InfN[FPDims[#],16]&/@inputrings;

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: FPDims of elements of fusion ring"
];

exportjson[ "fpdims", json ]


(* ::Subsection:: *)
(*FPdim (numerical test)*)


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = InfN[FPDim[#],16]&/@inputrings;

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: FPDim of fusion ring"
];

exportjson[ "fpdim", json ]


(* ::Subsection:: *)
(*group ring*)


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = GroupRingQ/@inputrings;

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: boolean stating whether ring is group ring"
];

exportjson[ "is_group_ring", json ]


(* ::Subsection:: *)
(*tensor product decompositions*)


inputrings = FRBC/@codes;
input = FC /@ inputrings;
output = WhichDecompositions/@inputrings/.r_FusionRing:>FC[r];

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: lists of formal codes of fusion rings arising in tensor product decompositions of fusion ring"
];

exportjson[ "tensor_product_decompositions", json ]


(* ::Subsection:: *)
(*adjoint fusion rings*)


inputrings = FRBC/@codes;
input = codes;
output = AdjointFusionRing/@inputrings/.r_FusionRing:>FC[r];

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: couple [ els, fc ] where els are the elements of the fusion ring that form the adjoint ring and fc is the formal code of the adjoint fusion ring"
];

exportjson[ "adjoint_fusion_rings", json ]


(* ::Subsection:: *)
(*upper central series*)


inputrings = FRBC/@codes;
input = codes;
output = UpperCentralSeries/@inputrings/.r_FusionRing:>FC[r];

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: upper central serie of fusion ring as list of couples [ [ els_1, fc_1 ], ..., [els_n, fc_n ] where els_i are the elements of the fusion ring with code fc_i that form the adjoint ring of fusion ring with code fc_{i-1}"
];

exportjson[ "upper_central_series", json ]


(* ::Subsection:: *)
(*adjoint irreps*)


inputrings = FRBC/@codes;
input = codes;
output = Sort[Sort/@AdjointIrreps[#]]&/@inputrings/.r_FusionRing:>FC[r];

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: partition of basis of fusion ring into adjoint irreps. Each subset is sorted and the partition is sorted lexicographically on the subsets"
];

exportjson[ "adjoint_irreps", json ]


(* ::Subsection:: *)
(*universal grading*)


inputrings = FRBC/@codes;
input = codes;
output = PMap[ UniversalGrading, inputrings]/.r_FusionRing:>FC[r];

json = tojson[ 
	input,
	output,
	"Input: formal code of fusion ring, " <> 
	"Output: partition of basis of fusion ring into adjoint irreps. Each subset is sorted and the partition is sorted lexicographically on the subsets"
];

exportjson[ "adjoint_irreps", json ]
