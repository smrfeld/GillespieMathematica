(* ::Package:: *)

BeginPackage["Gillespie`"]


fMakeReaction::usage="fMakeReaction makes a reaction dictionary, either unimolecular or 
bimolecular. Inputs in order: name, rate, reactants, products. Products is a list of
strings. For unimolecular reactions, reactant is a string; for bimolecular, reactants is a
list of strings."


fGillespie::usage="fGillespie runs the Gillespie algorithm. Inputs in order: tMax = max time.
tStEvery = timestep to write counts. cDict0 = association of initial particle counts.
uniRxnList = list of unimolecular reactions. biRxnList = list of bimolecular reactions.
consLaws = association of strings to counts of species to conserve. Returns {times, 
association of particle counts, association of reaction counts}."


Begin["Private`"]


fMakeReaction[name_,kr_,r_,p_]:=Return[<|"name"->name,"kr"->kr,"r"->r,"p"->p|>];


fGillespie[tMax_,tStEvery_,cDict0_,uniRxnList_,biRxnList_,consLaws_]:=Module[
{t,cDictSt,cDict,tSt,props,rxns,propsCum,prop,rxnChoose,dt,rxnListC,rxnNew,rxnCount,rxnCountSt}
,
(* Store counts *)
cDictSt=Association[];
Do[
If[!KeyExistsQ[consLaws,k],
cDictSt[k]={cDict0[k]};
];
,{k,Keys[cDict0]}];

(* Current counts *)
cDict=cDict0;

(* Write conservation for reagants into the reactions *)
rxnListC={};
Do[
rxnNew=rxn;
rxnNew["rC"]={};
rxnNew["r"]={};
Do[
	If[KeyExistsQ[consLaws,r],
		AppendTo[rxnNew["rC"],r],
		AppendTo[rxnNew["r"],r]
	];
,{r,rxn["r"]}];

(* Just throw out conserved products *)
If[rxn["p"]!={},
	rxnNew["p"]={};
	Do[
		If[!KeyExistsQ[consLaws,p],
			AppendTo[rxnNew["p"],p]
		];
	,{p,rxn["p"]}];
];
AppendTo[rxnListC,rxnNew];
,{rxn,Join[uniRxnList,biRxnList]}];

(* Store reaction counts *)
rxnCount=Association[];
rxnCountSt=Association[];
Do[
rxnCount[rxn["name"]]=0;
rxnCountSt[rxn["name"]]={0};
,{rxn,rxnListC}];

	(* Go over all time *)
t=0.0;
tSt={t};
Monitor[
While[t<tMax,

(* Propensities *)
	props={};
	rxns={};
	propsCum=0.0;

	(* Go through all possible reactions *)
	Do[

		(* Propensity *)
prop=Product[cDict[r],{r,rxn["r"]}];
prop*=Product[consLaws[r],{r,rxn["rC"]}];
prop*=rxn["kr"];

		(* Propensities *)
		propsCum+=prop;
		AppendTo[props,prop];
		AppendTo[rxns,rxn];
	,{rxn,rxnListC}];

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

(* Check if we need to store *)
If[IntegerPart[(t+dt)/tStEvery]==IntegerPart[t/tStEvery]+1,
(* Store *)
Do[
AppendTo[cDictSt[c],cDict[c]];
,{c,Keys[cDictSt]}];
AppendTo[tSt,(IntegerPart[t/tStEvery]+1)*tStEvery];
Do[
AppendTo[rxnCountSt[rxn["name"]],rxnCount[rxn["name"]]];
,{rxn,rxnListC}];
];

(* Advance time *)
t+=dt; 

(* Do the reaction *)
Do[
cDict[rs]-=1;
,{rs,rxnChoose["r"]}];
Do[
cDict[ps]+=1;
,{ps,rxnChoose["p"]}];

(* Count it *)
rxnCount[rxnChoose["name"]]+=1;
];
,
ProgressIndicator[t,{0,tMax}]
];

Return[{tSt,cDictSt,rxnCountSt}];
]


End[]


EndPackage[]
