(* ::Package:: *)

BeginPackage["PerformanceMetricFunctions`"]


(* The getMinDistance function calculates the minimum distance between inputConstellation1 symbols *)
getMinDistance[inputConstellation1_]:=Module[{DimMM1 = Dimensions[inputConstellation1]},
DMM = DistanceMatrix[DeleteDuplicates[ArrayReshape[inputConstellation1,{DimMM1[[1]]*DimMM1[[2]],DimMM1[[3]]}]]];
dmin = Min[Table[RankedMin[DMM[[xx,All]],2],{xx,1,Length[DMM]}]];
dmin
];

(* The getModulationQuantizer function calculates the CARDINALITY of the quantizer (non repeated values of DAC voltages) *)
getModulationQuantizer[inputConstellation2_]:=Module[{DimMM2 = Dimensions[inputConstellation2]},
DMM = DeleteDuplicates[ArrayReshape[inputConstellation2,{DimMM2[[1]]*DimMM2[[2]],DimMM2[[3]]}]];
QAM1 = Flatten[DMM[[All,{1}]]];
QAM2 = Flatten[DMM[[All,{2}]]];
QAM3 = Flatten[DMM[[All,{3}]]];
QAM4 = Flatten[DMM[[All,{4}]]];
QUANTIZER = Sort[DeleteDuplicates[Round[Join[QAM1,QAM2,QAM3,QAM4],0.0001]]];
QUANTIZER
];

(* The getModulationPerformanceMetrics function calculates performance metrics for the inputConstellation *)
getModulationPerformanceMetrics[inputConstellation_]:=Module[{DimMM = Dimensions[inputConstellation]},
(*Print["
\n
**********************************************************
 Starting the calculation of the performance metrics !!!
**********************************************************
\n
"];*)
(* OBS: The partitioned 4D constellation map to repeated symbols on the QAM partitions so the M = Length[ModulationMatrixSerial] calculation need to be considered only after the application of the "DeleteDuplicates" method. *)
ModulationMatrixSerial = DeleteDuplicates[ArrayReshape[inputConstellation,{DimMM[[1]]*DimMM[[2]],DimMM[[3]]}]];
M  = Length[ModulationMatrixSerial]; (* M = Number of Symbols *)
P = Table[Norm[ModulationMatrixSerial[[xx]]],{xx,1,M}]; (* P = summation of vectors *)
NumOfDim = DimMM[[3]];  (* Ndim = Dimension of modulation *)
DMM = DistanceMatrix[ModulationMatrixSerial];  (* DMM = distance matrix *)
d = Min[Table[RankedMin[DMM[[yy,All]],2],{yy,1,Length[DMM]}]];  (* d = minimum distance os symbols *)
Eave = Round[1/M Total[P],0.001];  (* Eave = Average symbol energy *)
Eb = Eave/Log2[M]; (* Eb = ENERGY per bit calculation *)
PEff = Round[d^2/(4Eb),0.001]; (* PEff = POWER efficiency calculation *)
SEff= Round[Log2[M]/(NumOfDim/2),0.001];  (* SEff = SPECTRAL efficiency calculation *)
cfm= Round[(d^2 NumOfDim)/(2*Eave),0.001]; (* cfm = CONSTELLATION FIGURE OF MERIT calculation *)
{distmin->d,AVE-> Eave, SE-> SEff,CFM-> cfm, CFMdB-> 10Log[10,cfm], PE-> PEff, PEdB-> 10Log[10,PEff], Dim-> NumOfDim, Number Of Symbols-> M}
];

plotQAMPartitions[QAMPartitions_,tokenprint_]:=Module[{pointsize = 0.04},
LEFTQAM = QAMPartitions[[1,All]];
RIGHTQAM = QAMPartitions[[2,All]];
gmpmL = getModulationPerformanceMetrics[LEFTQAM];
gmpmR = getModulationPerformanceMetrics[RIGHTQAM];

If[tokenprint==True,
dlq = Dimensions[LEFTQAM];
LEFTQAMRESHAPE=ArrayReshape[LEFTQAM,{dlq[[1]]*dlq[[2]],2}];
drq = Dimensions[RIGHTQAM];
RIGHTQAMRESHAPE=ArrayReshape[RIGHTQAM,{drq[[1]]*drq[[2]],2}];

tlq = Counts[LEFTQAMRESHAPE];
KEYLEFT = Keys[tlq];
VALLEFT = Values[tlq];
normalizationRatioL= 1/Max[VALLEFT];

LTLQ = Length[tlq];
ptlq = Table[Flatten[{KEYLEFT[[ii]],VALLEFT[[ii]]}],{ii,1,LTLQ}];
ptlq2 = Table[Tube[{Flatten[{KEYLEFT[[ii]],0}],Flatten[{KEYLEFT[[ii]],VALLEFT[[ii]]}]},0.05,VertexColors->{Green,Orange}],{ii,1,LTLQ}];
(*ptlq3 = Table[{Opacity[0.9],EdgeForm[None],Cylinder[{Flatten[{KEYLEFT[[ii]],0}],Flatten[{KEYLEFT[[ii]],VALLEFT[[ii]]}]},0.05]},{ii,1,LTLQ}];*)
functionPointsL = ListPlot3D[ptlq,ColorFunction->"Rainbow",PlotStyle->PointSize[0.04]];
pointsL = ListPointPlot3D[ptlq,ColorFunction->"Rainbow",PlotStyle->PointSize[0.04]];
(* PLOT OF LEFT QAM *)
Print[
{
Print[StringForm["\n********************* PLOT FOR LEFT QAM PARTITION ******************************\n"]],
Show[functionPointsL,ImageSize-> 350,PlotRange-> All,Axes-> True,Boxed-> True],
Show[pointsL,ImageSize-> 380,PlotRange-> All,Axes-> True,Boxed-> True],
Show[Graphics3D[ptlq2],ImageSize-> 480,PlotRange-> All,Axes-> True,Boxed-> True,BoxRatios-> {1,1,normalizationRatioL}]
(*,Show[Graphics3D[ptlq3],ImageSize\[Rule] 300,PlotRange\[Rule] All,Axes\[Rule] True,Boxed\[Rule] True,BoxRatios\[Rule] {1,1,normalizationRatioL},ViewPoint\[Rule] Default]*)
}
];

Print[gmpmL];

trq = Counts[RIGHTQAMRESHAPE];
KEYRIGHT = Keys[trq];
VALRIGHT = Values[trq];
normalizationRatioR = 1/Max[VALRIGHT];

LTRQ = Length[trq];
ptrq = Table[Flatten[{KEYRIGHT[[ii]],VALRIGHT[[ii]]}],{ii,1,LTRQ}];
ptrq2 = Table[Tube[{Flatten[{KEYRIGHT[[ii]],0}],Flatten[{KEYRIGHT[[ii]],VALRIGHT[[ii]]}]},0.05,VertexColors->{Green,Orange}],{ii,1,LTRQ}];
(*ptrq3 = Table[{Opacity[0.9],EdgeForm[None],Cylinder[{Flatten[{KEYRIGHT[[ii]],0}],Flatten[{KEYRIGHT[[ii]],VALRIGHT[[ii]]}]},0.05]},{ii,1,LTRQ}];*)
functionPointsR = ListPlot3D[ptrq,ColorFunction->"Rainbow",PlotStyle->PointSize[0.04]];
pointsR = ListPointPlot3D[ptrq,ColorFunction->"Rainbow",PlotStyle->PointSize[0.04]];
(* PLOT OF RIGHT QAM *)
Print[{
Print[StringForm["\n********************* PLOT FOR RIGHT QAM PARTITION ******************************\n"]],
Show[functionPointsR,ImageSize-> 350,PlotRange-> All,Axes-> True,Boxed-> True],
Show[pointsR,ImageSize-> 380,PlotRange-> All,Axes-> True,Boxed-> True],
Show[Graphics3D[ptrq2],ImageSize-> 480,PlotRange-> All,Axes-> True,Boxed-> True,BoxRatios-> {1,1,normalizationRatioR}]
(*,Show[Graphics3D[ptrq3],ImageSize\[Rule] 300,PlotRange\[Rule] All,Axes\[Rule] True,Boxed\[Rule] True,BoxRatios\[Rule] {1,1,normalizationRatioR},ViewPoint\[Rule] Default]*)
}
];

Print[gmpmR];
] (*end If*)
];

EndPackage[]






