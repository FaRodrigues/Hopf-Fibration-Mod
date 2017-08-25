(* ::Package:: *)

BeginPackage["HopFibrationFunctions`"]

ClearAll[\[Theta]1,\[Phi]1,\[Theta]x,\[Phi]x];
ClearAll[splineCircle,pointsAndConnection];

StokesToPolarCoordinates::usage = "Realize the conversion from normalized Stokes vectors angular coordinates";
getHopfFiber::usage = "Returns the Inverse Hopf map";
rightDIRECTHOPFMAP::usage = "Returns the Right Direct Hopf map";
factorProjSTEREO::usage = "Returns the distorted stereographic projection";
getRGBColors::usage = "Returns RGB colors based in parameters";

Begin["\.b4Private`"]
(* This part of code was obtained from Internet *)
(*------https://mathematica.stackexchange.com/questions/10957/an-efficient-circular-arc-primitive-for-graphics3d/10994 ---------*)

splineCircle[m_List,r_,angles_List: {0,2 \[Pi]}]:=Module[{seg,\[Phi],start,end,pts,w,k},{start,end}=Mod[angles//N,2 \[Pi]];
If[end<=start,end+=2 \[Pi]];
seg=Quotient[end-start//N,\[Pi]/2];
\[Phi]=Mod[end-start//N,\[Pi]/2];
If[seg==4,seg=3;\[Phi]=\[Pi]/2];
pts=r RotationMatrix[start].#&/@Join[Take[{{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1}},2 seg+1],RotationMatrix[seg \[Pi]/2].#&/@{{1,Tan[\[Phi]/2]},{Cos[\[Phi]],Sin[\[Phi]]}}];
If[Length[m]==2,pts=m+#&/@pts,pts=m+#&/@Transpose[Append[Transpose[pts],ConstantArray[0,Length[pts]]]]];
w=Join[Take[{1,1/Sqrt[2],1,1/Sqrt[2],1,1/Sqrt[2],1},2 seg+1],{Cos[\[Phi]/2],1}];
k=Join[{0,0,0},Riffle[#,#]&@Range[seg+1],{seg+1}];
BSplineCurve[pts,SplineDegree->2,SplineKnots->k,SplineWeights->w]]/;Length[m]==2||Length[m]==3

pointsAndConnection[points_]:=Sequence@@{Sequence@@Point/@#,Line@#}&@points
surroundingCircles=GeometricTransformation[splineCircle[{0,0,0},1],{{RotationMatrix[0,{1,0,0}],{0,0,0}},{RotationMatrix[Pi/2,{1,0,0}],{0,0,0}},{RotationMatrix[Pi/2,{0,1,0}],{0,0,0}}}];
poincareSphereList={{White,Opacity@0.3,Sphere[{0,0,0},1],Opacity@1,Thickness@0.004,PointSize@0.01,Black,pointsAndConnection@{{0,0,1},{0,0,-1}},Black,pointsAndConnection@{{1,0,0},{-1,0,0}},Black,pointsAndConnection@{{0,1,0},{0,-1,0}},Black,Point[{0,0,0}],Black,Thin,surroundingCircles},Boxed->False,PlotRange->ConstantArray[{-1,1},3],ImageSize->300,RotationAction->"Clip"};

(* End of code obtained from Internet *)
(*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------*)
(* 
StokesToPolarCoordinates function is an implementation of the conversion between Stokes and angular parameters (spherical coordinates).
StokesVector_ input is the normalized Stokes vector representation with (s1,s2,s3) parameters;
column_ input is used to select whether the function returns Stokes or angular parameters
The Stokes parameters output can be used, for example, to calculate the amount of ERROR inserted betwwen vector entries and outputs
*)
StokesToPolarCoordinates[StokesVector_,column_]:= Module[
{radius = Sqrt[(StokesVector[[1]])^2+(StokesVector[[2]])^2+(StokesVector[[3]])^2]},
ZDFlag = 10^-16; (* ZDFlag is used to avoid division by zero - CAUTION this cumbersome trick insert ERROR on the calculation *)
RS0 = Solve[N[ArcCos[StokesVector[[3]]/(radius+ZDFlag)]]==\[Theta]1, \[Theta]1,Reals] //Flatten; (* Angle \[Theta]1 determines latitude in the S2 sphere *)
RS4 = Solve[N[ArcTan[StokesVector[[2]]/(StokesVector[[1]]+ZDFlag)]]==\[Phi]1, \[Phi]1,Reals] //Flatten; (* Angle \[Phi]1 determines longitude in the S2 sphere *)
RS5 = Flatten[Join[RS0,RS4]];
If[
(StokesVector[[1]]) <0, 
(* Specific rule for \[Theta]x calculation *)
	bas1 = {\[Phi]x->\[Phi]1, \[Theta]x->-\[Theta]1}  /.RS5,
	bas1 = {\[Phi]x->\[Phi]1, \[Theta]x->\[Theta]1}  /.RS5
];
outstokes = Round[{radius*Sin[\[Theta]x]Cos[\[Phi]x],radius*Sin[\[Theta]x]Sin[\[Phi]x],radius*Cos[\[Theta]x]} /.bas1,0.0001];
outangular = {\[Phi]->\[Phi]x,\[Theta]->\[Theta]x} /.bas1;
out={outstokes,outangular};
out[[column]]  //Simplify
];

(* 
getHopfFiber function is an example of implementation of the "Inverse Hopf Map Equation" .
PHI,TH and \[Psi]h are respectively the input parameters for \[Phi], \[Theta] and \[Psi] angles
*)
getHopfFiber[{PHI_,TH_,PSI_}, SA_]:=Module[{psh=PSI},
If[And[Abs[PHI]==0,If[Abs[TH]!=0,TH-1.5707<0.0001]],psh=PSI-\[Pi]/4,psh=PSI];
z1 = Cos[TH/2]Exp[I((psh+SA)+PHI/2)];
z2 = Sin[TH/2]Exp[I((psh+SA)-PHI/2)];
Cz1 = Simplify[ComplexExpand[ReIm[z1]]];
Cz2 = Simplify[ComplexExpand[ReIm[z2]]];
Round[Flatten[{Cz1,Cz2}],0.0001]
];

(* 
rightDIRECTHOPFMAP function is an implementation of the "Direct Hopf Map Equation" .
{x_,y_,z_,w_} are the coordinates of the 4D vector (Quaternion) input
{S1,S2,S3} are the coordinates of the Stokes vector output
*)
rightDIRECTHOPFMAP[{x_,y_,z_,w_}]:=Module[{S0=x^2+y^2+z^2+w^2},
S1 = x^2+y^2-z^2-w^2;
S2 = 2*((y*w)+(x*z));
S3 = 2*((x*w)-(y*z));
If[S0 <= 10^-9,
stokes = {S1,S2,S3},
stokes = {S1,S2,S3}/S0
];
stokes
];

(* 
factorProjSTEREO function is an implementation of the "Stereographic Projection" .
{x_,y_,z_,w_} are the coordinates of the 4D vector (Quaternion) input
factor_ is used to distort the standard stereographic projection (used only for vizualization purposes)
r is used to avoid division by zero
*)

factorProjSTEREO[{x_,y_,z_,w_},factor_] := Module[{r=N[1-factor]},
X=x/(1-(r*w));Y=y/(1-(r*w));Z=z/(1-(r*w));
{X,Y,Z}
];

getRGBColors[{a_,b_,c_}]:=Module[{A=N[(1+a)/2],B=N[(1+b)/2],C=N[(1+c)/2]},
RGBColor[A,B,C]
];

basisPoints = {};

setBasisPointsFromTXTFile[namefile_]:= Module[{filepath = FileNameJoin[{NotebookDirectory[],"basispoints",namefile}]},
Print[StringForm["\nLoading basis points from '``'.\n",filepath]];
 icosphere42M90=Import[filepath,"Table"];
 basisPoints = DeleteDuplicates[Round[icosphere42M90,0.001]]
 ]

get4DConstellationFromSampledHopfFibration[bpoints_,PSIsamples_,startAngle_]:=Module[{bp=bpoints},
If[Length[bp]>0,
Print[StringForm["\nStokes vectors ENTRIES:\n``",bp]];
StokesCoordinatesOutput =Table[StokesToPolarCoordinates[bp[[scoindex]],1],{scoindex,1,Length[bp]}];
Print[StringForm["\nStokes vectors OUTPUTS:\n``",StokesCoordinatesOutput]];
(* 
Calculates the Euclidean distance between basePoints and StokesCoordinatesOutput. 
This distance computes the ERROR inserted by the conversion StokesToPolarCoordinates.
*)
DIST =Round[ EuclideanDistance[bp,StokesCoordinatesOutput],0.0000001];
Print[StringForm["\nEuclidean distance between Stokes vectors ENTRIES and OUTPUTS:\n``",DIST]];
PolarCoordinatesOutput =Table[StokesToPolarCoordinates[bp[[pcoindex]],2],{pcoindex,1,Length[bp]}];
Print[StringForm["\nStokes vectors in the ANGULAR parametric representation:\n``",PolarCoordinatesOutput]];


ModulationMatrix = Table[
Table[discreteparam={\[Phi],\[Theta],\[Psi]}/.PolarCoordinatesOutput[[hh]];getHopfFiber[discreteparam,startAngle],
{\[Psi],Range[0,2\[Pi](1-1/PSIsamples),(2\[Pi])/PSIsamples]}],
{hh,1,Length[bp]}];

DimM = Dimensions[ModulationMatrix];
NumberOfConstellationSymbols = DimM[[1]]*DimM[[2]];

Print[StringForm["\nThe Four-Dimensional constellation has `` symbols, obtained from 'Sampled Discrete Hopf Fibration' with `` base points and `` angular samples.",NumberOfConstellationSymbols,DimM[[1]],DimM[[2]]]];

(* The numberOfPlotListPoints variable defines the number of points in each Hopf circle *)
numberOfPlotListPoints = 360; 

ModulationMatrixContinuous = Table[
Table[param={\[Phi],\[Theta],\[Psi]}/.PolarCoordinatesOutput[[hh]];getHopfFiber[param,startAngle],
{\[Psi],Range[0,2\[Pi](1-1/numberOfPlotListPoints),(2\[Pi])/numberOfPlotListPoints]}],
{hh,1,Length[bp]}];

{ModulationMatrix,ModulationMatrixContinuous} (* get4DConstellationFromSampledHopfFibration OutPut *)
,
Print["ERROR MESSAGE:\nBasis points are NOT defined !!! Please call setBasisPointsFromTXTFile before calling get4DConstellationFromSampledHopfFibration"];
]
];
(* End of get4DConstellationFromSampledHopfFibration *)
(*-------------------------------------------------------------------------------------*)

getQAMFrom4DModulation[ModulaMatrix_]:=Module[{i1=1,i2=2,i3=3,i4=4},
LEFTQAM = ModulaMatrix[[All,All,i1;;i2]]; (* Left 2D (QAM) Partition *)
RIGHTQAM = ModulaMatrix[[All,All,i3;;i4]];     (* Right 2D (QAM) Partition *)
NDLEFTQAM = DeleteDuplicates[LEFTQAM];
NDRIGHTQAM=DeleteDuplicates[RIGHTQAM];
Print[
{
Show[ListPlot[NDLEFTQAM,PlotStyle-> PointSize[0.04],ImageSize-> 250,AspectRatio->Automatic, PlotLabel->"LEFT QAM"]],
Show[ListPlot[NDRIGHTQAM,PlotStyle-> PointSize[0.04],ImageSize-> 250,AspectRatio->Automatic, PlotLabel->"RIGHT QAM"]],
Show[ListPlot[Join[NDLEFTQAM,NDRIGHTQAM],PlotStyle-> PointSize[0.04],ImageSize-> 250,AspectRatio->Automatic, PlotLabel->"SUPERPOSED QAM"]]
}
];
{LEFTQAM,RIGHTQAM}
];

(*-------------------------------------------------------------------------------------*)

plotHopfFibration[basePointsInput_,ModulationMatrixInput_,ModulationMatrixContinuousInput_,samples_,FactorProj_]:=Module[{diMMC = Dimensions[ModulationMatrixContinuous]},
ModulationMatrixHOPFPoints = Table[
Table[{getRGBColors[basePointsInput[[mm1]]],PointSize[0.04],Point[rightDIRECTHOPFMAP[ModulationMatrixInput[[mm1,vv1]]]]},{vv1,1,samples}],
{mm1,1,Length[ModulationMatrix]}
];

ModulationMatrixSTEREOPoints = Table[
Table[{getRGBColors[basePointsInput[[mm2]]],PointSize[0.04],Point[factorProjSTEREO[ModulationMatrixInput[[mm2,vv2]],FactorProj]]},{vv2,1,samples}],
{mm2,1,Length[ModulationMatrix]}
];

ModulationMatrixSTEREOContinuous = Table[
Table[{getRGBColors[basePointsInput[[mm3]]],PointSize[0.015],Point[factorProjSTEREO[ModulationMatrixContinuousInput[[mm3,vv3]],FactorProj]]},{vv3,1,diMMC[[2]]}],
{mm3,1,diMMC[[1]]}
];
Print[
{
Show[Graphics3D[Join[{ModulationMatrixSTEREOPoints,ModulationMatrixSTEREOContinuous,Opacity[.01],Sphere[{0,0,0}]}], PlotLabel->"Stereographic projection - Inverse Hopf Map"],ImageSize-> 400,PlotRange-> All,Axes-> False,Boxed-> False,ViewPoint-> {Pi,0,0.2}],
Show[Graphics3D[Join[{ModulationMatrixSTEREOPoints,ModulationMatrixSTEREOContinuous,Opacity[.01],Sphere[{0,0,0}]}], PlotLabel->"Stereographic projection - Inverse Hopf Map"],ImageSize-> 350,PlotRange-> All,Axes-> False,Boxed-> True,ViewPoint-> Top],
Show[Graphics3D[Join[{ModulationMatrixHOPFPoints,poincareSphereList,Opacity[.1],Sphere[{0,0,0}]}], PlotLabel->"Direct Hopf Map"],ImageSize-> 400,PlotRange-> All,Axes-> False,Boxed-> False,ViewPoint-> {Pi,0.35,0.35}]
}
];
];
End[]
EndPackage[];



