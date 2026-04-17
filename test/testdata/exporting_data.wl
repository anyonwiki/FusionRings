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
	Export["/home/gert/Projects/FusionRings/test/data/"<>fn<>".json",dt,"JSON"];


(* ::Section:: *)
(*Creation*)


(* ::Text:: *)
(*Check that creation functions of FusionRings.jl give same result as those of Anyonica*)


(* ::Input:: *)
(*zntables = *)
(*tojson[*)
(*ToString/@Range[8],*)
(*MT[FusionRingZn[#]]&/@Range[8],*)
(*"Input: int i, Output: multiplication table for ring of cyclotomic group of order i"*)
(*];*)
(*exportjson["zn_tables",zntables]*)


(* ::Input:: *)
(*smallgroups = *)
(*Join[Flatten[*)
(*Table[ {"AbelianGroup",{i,j}},{i,2,3},{j,2,3} ],1],*)
(*Table[{"DihedralGroup",i}, {i,3,5 }],*)
(*Table[{"SymmetricGroup",i}, {i,4 }]*)
(*];*)
(*smallgroupsmulttabs = FiniteGroupData[#,"MultiplicationTable"]&/@smallgroups;*)
(**)
(*grouptables = *)
(*tojson[*)
(*smallgroupsmulttabs,*)
(*MT@*FusionRingFromGroup/@smallgroupsmulttabs,*)
(*"Input: mult tab of small group, Output: mult tab of fusion ring of that group"*)
(*];*)
(*exportjson["grouptables",grouptables];*)


(* ::Input:: *)
(*smallabeliangroups = Flatten[Table[ FiniteGroupData[{"AbelianGroup",{i,j}},"MultiplicationTable"],{i,2,3},{j,2,3} ],1];*)
(*hitables = *)
(*	tojson[*)
(*	smallabeliangroups,*)
(*	MT@*FusionRingHI/@smallabeliangroups,*)
(*"Input: mult tab of small abelian group, Output: mult tab of HI fusion ring of that group"*)
(*];*)
(*exportjson["HI_tables",hitables]*)


(* ::Input:: *)
(*tytables =*)
(*tojson[*)
(*smallgroupsmulttabs,*)
(*MT@*FusionRingTY/@smallgroupsmulttabs,*)
(*"Input: mult tab of small group, Output: mult tab of TY fusion ring of that group"*)
(*];*)
(*exportjson["TY_tables",tytables]*)


(* ::Input:: *)
(*psu2ktables = *)
(*tojson[*)
(*Range[8],*)
(*MT@*FusionRingPSU2k /@ Range[8],*)
(*"Input: int i, Output: mult tab of PSU(2)_i fusion ring"*)
(*];*)
(**)
