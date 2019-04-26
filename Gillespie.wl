(* ::Package:: *)

BeginPackage["Gillespie`"]


makeReaction::usage="makeReaction makes a reaction dictionary, either unimolecular or 
bimolecular. Inputs in order: name, rate, reactants, products. Reactions and products are lists of
strings."


runGillespie::usage="runGillespie runs the Gillespie algorithm. Inputs in order: 
tMax = max time. 
countsInit = association of initial particle counts.
rxnList = list of reactions. 
conservedSpecies = list of species to conserve, else empty. 
drWrite = If specified, write out at this timestep; if None, write only the times when reactions occur 
or particle counts change
Returns:
association of species -> {times,particle counts}
association of reactions -> {times,reaction counts}."


Begin["Private`"]


makeReaction[name_,kr_,r_,p_]:=Module[
{rxn},
	If[Length[r]==0&&Length[p]==0,
		Print["Empty reactions are forbidden"];
		Return[None];
	];
	If[Length[r]==0&&Length[p]!=1,
		Print["Only 0->A is allowed, no further products."];
		Return[None];
	];
	rxn=<|"name"->name,"kr"->kr,"r"->r,"p"->p|>;
	rxn["species"]=DeleteDuplicates[Join[r,p]];
	Return[rxn];
];


runGillespie[tMax_,countsInit_,allRxnsList_,conservedSpecies_,dtWrite_:None]:=Module[
{t,countsSt,countsCurrent,props,rxns,propsCum,prop,rxnChoose,dt,rxnNew,rxnCountSt,iStep,iStepNew,rxnCountsCurrent}
,
(* Store counts *)
countsSt=Association[];
Do[
	countsSt[k]={{0.0,countsInit[k]}};
,{k,Keys[countsInit]}];

(* Current counts *)
countsCurrent = countsInit;

(* Store reaction counts *)
rxnCountSt=Association[];
rxnCountsCurrent=Association[];
Do[
	rxnCountSt[rxn["name"]]={{0.0,0}};
	rxnCountsCurrent[rxn["name"]]=0;
,{rxn,allRxnsList}];

(* Go over all time *)
t=0.0;
iStep=0;
Monitor[
	While[t<tMax,

		(* Propensities *)
		props={};
		rxns={};
		propsCum=0.0;

		(* Go through all possible reactions *)
		Do[

			(* Propensity *)
			(* Note: if type 0\[Rule]A, make prop to A *)
			If[Length[rxn["r"]]!=0,
				prop=Product[countsCurrent[r],{r,rxn["r"]}],
				prop=countsCurrent[rxn["p"][[1]]]
			];
			prop*=rxn["kr"];

			(* Propensities *)
			propsCum+=prop;
			AppendTo[props,prop];
			AppendTo[rxns,rxn];
			
		,{rxn,allRxnsList}];

		If[Total[props]==0.0,
			(* No reactions left; finish *)
			Break[];
		];

		(* Choose a reaction *)
		rxnChoose=RandomChoice[props->rxns];
		
		(* Time of next reaction *)
		dt=Log[1.0/RandomReal[{0,1}]]/propsCum;

		(* Check this time is still ok *)
		If[t+dt>tMax,Break[]];

		(* Advance time *)
		t+=dt; 
		If[dtWrite=!=None,
			iStepNew = IntegerPart[t/dtWrite];
			If[iStepNew!=iStep,
				(* Write *)
				Do[
					Do[
						AppendTo[countsSt[species],{iStepWrite*dtWrite,countsCurrent[species]}];
					,{species,Keys[countsSt]}];
					Do[
						AppendTo[rxnCountSt[rxn],{iStepWrite*dtWrite,rxnCountsCurrent[rxn]}];
					,{rxn,Keys[rxnCountsCurrent]}];
				,{iStepWrite,iStep+1,iStepNew}];
				(* Advance *)
				iStep=iStepNew;
			];
		];

		(* Do the reaction *)
		rxnCountsCurrent[rxnChoose["name"]]+=1;
		Do[
			If[!MemberQ[conservedSpecies,rs],
				countsCurrent[rs]-=1;
			];
		,{rs,rxnChoose["r"]}];
		Do[
			If[!MemberQ[conservedSpecies,ps],
				countsCurrent[ps]+=1;
			];
		,{ps,rxnChoose["p"]}];
		
		(* Store counts at the exact times *)
		If[dtWrite===None,
			Do[
				AppendTo[countsSt[species],{t,countsCurrent[species]}];
			,{species,rxnChoose["species"]}];
			AppendTo[rxnCountSt[rxnChoose["name"]],{t,rxnCountsCurrent[rxnChoose["name"]]}];
		];
	];
,ProgressIndicator[t,{0,tMax}]];

Return[{countsSt,rxnCountSt}];
]


End[]


EndPackage[]
