(* ::Package:: *)

(* ::Title:: *)
(*Test  Data  for  FusionRings  Package*)


<<Anyonica`
codes = Import["/home/gert/Projects/FusionRings/test/data/test_fusion_ring_codes.wdx"];


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

input = Transpose @ { MT@*FRBC/@randomcodes, perms }; 

output = 
	Table[ 
		MT[PermutedRing[FusionRing[ "MultiplicationTable"->inp[[1]]],#]]&/@inp[[2]],
		{inp,input}
	];

permtabs = 
	tojson[ 
		input,
		output,
		"Input: list { mt, perms } where mt is original mult tab of a ring and perms is list of 2 permutations, Output: list of mult tabs corresponding to permutations of mt"
	];
	
exportjson[ "permuted_tabs", permtabs ]


(* tensor products between rings. We can only test for small rings *)
maxrank4codes = Cases[ codes, c_ /; c[[1]] < 5 ];

{ multfreecodes, multcodes } = BinSplit[ maxrank4codes, #[[2]] == 1& ];

SeedRandom[1];

randomcodes = Join[ RandomSample[ multfreecodes, 4 ], RandomSample[ multcodes, 4 ] ];

input = 
	Map[ FRBC, Tuples[{ randomcodes, randomcodes }], {2} ];

output = PMap[ MT[TensorProduct@@#]&, input ];

tptabs = tojson[ 
	Map[ MT, input, {2} ],
	output,
	"Input: couple of mult tabs of 2 fusion rings, Output: mult tab of tensor product of the rings"
];

exportjson[ "tensor_product_tables", tptabs ]


inputrings = FRBC/@codes;
input = MT /@ inputrings;
output = SubFusionRings/@inputrings/.r_FusionRing:>MT[r];

subrings = tojson[ 
	input,
	output,
	"Input: mult tab of fusion ring, Output: either (1) empty list if ring has no non-trivial sub fusion rings or (2) list of couples ( els, mt ) where els are the elements in the parent ring that make up the subring and mt is the multiplication table of the subring"
];

exportjson[ "sub_ring_data", subrings ]


inputrings = FRBC/@codes;
input = MT /@ inputrings;
output = Sort@*Rest@*FRA/@inputrings;

autos = tojson[ 
	input,
	output,
	"Input: mult tab of fusion ring, Output: non-trivial fusion ring automorphisms as permutation vectors, sorted lexicographically"	
];

exportjson[ "automorphisms", autos ]
