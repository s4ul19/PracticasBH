(* ::Package:: *)

BeginPackage["MatrixComponent`"]
Basis::usage="Basis[n_,m_]\nGenera la Base de Fock de una configuracion unidimensional de n=NumeroDeParticulas y m=NumeroDeSitios sitios";
hamiltonianBH::usage="hamiltonianBH[NumeroDeParticulas_,NumeroDeSitios_,ConstanteCinetica_,ConstanteInteraccion_]\n Crea el Hamiltoniano de BoseHubbard unidimensional, segun una l\[OAcute]gica de primeros vecinos.";
DeepHamiltonianBH::usage="M\[AAcute]s optimo que hamiltonianBH\nDeepHamiltonianBH[NumeroDeParticulas_,NumeroDeSitios_,ConstanteCinetica_,ConstanteInteraccion_]\n Crea el Hamiltoniano de BoseHubbard unidimensional, segun una l\[OAcute]gica de primeros vecinos"
CreateArray::usage="CreateArray[position_,values_]\nCrea Sparsearray segun una lista de posiciones y una lista de valores.\nLos objetos de la lista de posiciones deben ser de la forma {i,j}, i y j enteros\n position[[i]] le correspende values[[i]]. "
EigenInfo::usage="EigenInfo[sparse_]\n Se le pasa un SparseArray y retorna la informacion sobre sus Eigenvectores (Primera entrada) y sobre sus respectivos eigenvalores (Segunda entrada)."
UnitarySparse::usage="UnitarySparse[OldBase_,NewBase_]\nDada dos bases de un mismo espacio vectorial, retorna la matriz de cambio de base."
SplitSymmetricBasis::usage="SplitSymmetricBasis[basis_]\n Retorna un Diccionario de la base sim\[EAcute]trica asociada a la base dada.\n Mediante SplitSymmetricBasis[''Symmetric''] se obtiene los elementos sim\[EAcute]tricos y con SplitSymmetricBasis[''Antisymmetric''] los elementos antisim\[EAcute]tricos."
ParityRepresentationHamiltonian::usage="ParityRepresentationHamiltonian[n_,m_,J_,U_]\n Retorna la representacion sectorizada bajo la simetr\[IAcute]a de paridad del hamiltoniano de BoseHubbard."
SymmetricSectorHamiltonian::usage="SymmetricSectorHamiltonian[n_,m_,J_,U_,sector_]\nRetorna el sector sim\[EAcute]trico con sector=''Symmetric''\n Retorna el sector antisim\[EAcute]trico con sector=''Antisymmetric''\n Hamiltoniano de Bose Hubbard "
KLevelSpacing::usage="KLevelSpacing[list_,order_]"
DataExtraction::usage="DataExtraction[route_,type_]\n Extrae las archivos tipo type de un directorio localizado en route."
CreateHistogram::usage="CreateHistogram[list_]\nCrea histograma PDF seg\[UAcute]n una lista de datos."
DiscretePDF::usage="DiscretePDF[list_]\nRetorna la distribuci\[OAcute]n discreta de probabilidad, Retorna duplas, primera entrada coordenada x, representa valor representativo del bin ; segunda coordenada y, representa la probabilidad asociada al bin"
Spacing::usage="Spacing[list_,order_] \n Crea una lista de los espaciamientos relativos de orden n de la lista de datos."
Begin["`Private`"]
DiscretePDF[list_]:=Module[{L=HistogramList[list,"FreedmanDiaconis","PDF"]},
Table[{L[[1,i]]+(L[[1,i+1]]-L[[1,i]])/2,L[[2,i]]},{i,Length[L[[2]]]}]]   

Spacing[list_,order_]:=Module[{},Table[{list[[i+order]]-list[[i]]},{i,1,Length[list]-order}]]

DataExtraction[route_,type_]:=Module[{files=FileNames["*."<>type,route],DataList},
DataList=Import/@files;
AssociationThread[FileNameTake/@files->DataList]]

listaBase[v_,n_]:=Table[If[i==1,v,0],{i,n}]; 

Indice[lista_]:=Module[{n=Length[lista]},SelectFirst[Range[n-1],AllTrue[lista[[#+1;;n-1]],(#==0)&]&,Missing["NotFound"]]]

BaseN[lista_,n_]:=Module[{k=Indice[lista],copia=lista},
copia[[k]]-=1;
If[k+1<=Length[copia],copia[[k+1]]=n-Total[copia[[1;;k]]];];
Do[copia[[i]]=0,{i,k+2,Length[copia]}];
copia]

Repetir[f_,ini_,n_]:=NestList[f,ini,n-1]
Basis[N_,M_]:=Module[{Base1=listaBase[N,M],n1=(M+N-1)!/(N!*(M-1)!)},Repetir[BaseN[#,N]&,Base1,n1]]

Tag[lista_]:= Module[{n=Length[lista]},Sum[Sqrt[100.*i+3]*lista[[i]],{i,n}]] 


(*Calcula el arreglo de tags*)
ArrayTag[lista_]:=Module[{n=Length[lista],A1={}},
Do[AppendTo[A1,Tag[lista[[i]]]],{i,n}];A1]

Posicion[valor_,lista_]:=Module[{pos=Position[lista,valor]},If[pos==={},Return[Null],First[pos]]](*Detecta si el valor dado esta en la lista y si lo esta da la posicion en la que se halla.*)




singleMatrizAij[listB_,listTag_]:=Module[{v1=Length[listB],A=AplicarOperatorAij[listB],k1=Posicion[Tag[listB],listTag]},
Table[{{Posicion[Tag[A[[i]][[1]]],listTag],k1},A[[i]][[2]]},{i,Length[A]}]]


CrearSparseArrayAij[lista_,index_]:=Module[{pares},pares={Flatten[#[[1]]],#[[2]]}&/@lista;
SparseArray[Rule@@@pares,index]]


TotalSparseAij[listaB_,listaTag_]:=Sum[CrearSparseArrayAij[singleMatrizAij[listaB[[i]],listaTag],Length[listaB]],{i,1,Length[listaB]}]




OperatorNi[lista_,entrada_]:=Module[{k=entrada,copy=lista},
{copy,copy[[k]]*(copy[[k]]-1)}]



AplicarOperatorNi[lista_]:=Module[{resultados,copia,suma},resultados=Table[OperatorNi[lista,k],{k,Length[lista]}];
copia=resultados[[1,1]];suma=Total[resultados[[All,2]]];{copia,suma}]


singleMatrizNi[listB_,listTag_]:=Module[{v1=Length[listB],A=AplicarOperatorNi[listB],k1=Posicion[Tag[listB],listTag]},
{{k1,k1},A[[2]]}]


CrearSparseArrayNi[list_,index_]:=Module[{k=Flatten[list[[1]]],I=index,A=list[[2]]},SparseArray[{k->A},{I,I}]]


TotalSparseNi[listaB_,listaTag_]:=
Sum[CrearSparseArrayNi[singleMatrizNi[listaB[[i]],listaTag],Length[listaB]],{i,1,Length[listaB]}]


HKin[Sparse_,val_]:=Module[{J=val},
-(J/2)*(ConjugateTranspose[Sparse]+Sparse)]


HInt[Sparse_,val_] :=Module[{U=val},
(U/2)*Sparse]

hamiltonianBH[n_,m_,J_,U_]:=Module[{U1=U,J1=J},
V=Basis[n,m];
TagV=Sort[ArrayTag[Basis[n,m]]];
TKin=TotalSparseAij[V,TagV] ;
TInt=TotalSparseNi[V,TagV];
H1=HKin[TKin,J1];
H2=HInt[TInt,U1];
H1+H2 ]

OperatorAij[i_,j_,lista_]:=Module[{copy=lista,valorI,valorII},If[i==j,Return[Null]];
valorI=copy[[i]];
valorII=copy[[j]];
copy[[i]]=valorI+1;
copy[[j]]=valorII-1;
If[Min[copy]<0,Return[Null]];
{copy,Sqrt[(valorI+1.)*valorII]}]


AplicarOperatorAij[lista_]:=Module[{n=Length[lista],pares},pares=Flatten[Table[{{i,i-1},{i-1,i}},{i,2,n}],1];
DeleteCases[OperatorAij[#[[1]],#[[2]],lista]&/@pares,Null]]


DeepHamiltonianBH[n_,m_,J_,U_]:=Module[{V,assoc,len,b,neighbors,v,coef,j0,rulesKin,TKin,diag,H2},(*Genera base y asociaci\[OAcute]n vector->\[IAcute]ndice*)V=Basis[n,m];
len=Length[V];
assoc=AssociationThread[V->Range[len]];
rulesKin=Reap[Do[b=V[[i]];
neighbors=AplicarOperatorAij[b];
Do[{v,coef}=neighbors[[j]];
j0=assoc[v];
Sow[{i,j0}->coef],{j,1,Length[neighbors]}],{i,1,len}]][[2,1]];
TKin=If[rulesKin==={},SparseArray[{},{len,len}],SparseArray[rulesKin,{len,len}]];
diag=(U/2)*(Total[#*(#-1)]&/@V);


H2=SparseArray[Band[{1,1}]->diag,{len,len}];

-J*TKin+H2]



CreateArray[pos_,val_]:=Module[{Asoc=AssociationThread[pos,N[val]]//Normal,n=Max[Flatten[pos]]},SparseArray[Asoc,{n,n}]]

EigenInfo[sparse_]:=Module[{eigensys=Eigensystem[N[sparse]]},Transpose[{Normalize/@eigensys[[2]],eigensys[[1]]}]]

UnitarySparse[basis1_,basis2_]:=Module[{EA=basis1,EB=basis2,n=Length[basis1]},If[Length[EA]!=Length[EB],Message[UnitarySparse::bdim];
Return[$Failed];];
SparseArray[Outer[Dot,N[EB],N[EA],1]]]



SplitSymmetricBasis[list_]:=Module[{basis=list,symPairs,palindromes,symMap,antiMap,palMap},palindromes=Select[basis,#===Reverse[#]&];
palMap=Map[{#,#,1,0}&,palindromes];
symPairs=DeleteDuplicatesBy[Select[basis,#=!=Reverse[#]&],Sort[{#,Reverse[#]}]&];
symMap=Map[Function[{v},With[{w=Reverse[v]},If[OrderedQ[{v,w}],{v,w,1./Sqrt[2],1./Sqrt[2]},{w,v,1./Sqrt[2],1./Sqrt[2]}]]],symPairs];
antiMap=Map[Function[{v},With[{w=Reverse[v]},If[OrderedQ[{v,w}],{v,w,1./Sqrt[2],-1./Sqrt[2]},{w,v,1./Sqrt[2],-1./Sqrt[2]}]]],symPairs];
<|"Symmetric"->Join[palMap,symMap],"Antisymmetric"->antiMap|>]

ScaledCanonicalVectorij[i_,j_,vali_,valj_,dim_]:=Module[{id=ConstantArray[0,dim]},id[[i]]+=vali;
If[i!=j,id[[j]]+=valj];
id]

ParityRepresentationHamiltonian[n_,m_,J_,U_]:=Module[{H,B,len,assoc,SymmetricBasis,PartSymmetric,PartAntisymmetric,Entries,EntriesPartTotal,NewBasis,U1,mat,threshold,NewH1},H=DeepHamiltonianBH[n,m,J,U];
B=Basis[n,m];
len=Length[B];
SymmetricBasis=SplitSymmetricBasis[B];
PartSymmetric=SymmetricBasis["Symmetric"];
PartAntisymmetric=SymmetricBasis["Antisymmetric"];
assoc=AssociationThread[B,Range[len]];
Entries=Join[PartSymmetric,PartAntisymmetric];
EntriesPartTotal=Entries/. {a_,b_,x_,y_}:>{{assoc[a],assoc[b]},x,y};
NewBasis=Table[ScaledCanonicalVectorij[#1,#2,#3,#4,len]&@@{EntriesPartTotal[[i,1,1]],EntriesPartTotal[[i,1,2]],EntriesPartTotal[[i,2]],EntriesPartTotal[[i,3]]},{i,1,Length[EntriesPartTotal]}];
U1=SparseArray[NewBasis];
mat=U1 . H . Transpose[U1];
NewH1=Chop[mat,10^-10]  (*Usar Chop con tolerancia adecuada*)]

SymmetricSectorHamiltonian[n_,m_,J_,U_,sector_]:=Module[{Horig,B,len,assoc,SymmetricBasis,entries,newBasis,U1,Hsec},Horig=DeepHamiltonianBH[n,m,J,U];(*Hamiltoniano en base original*)B=Basis[n,m];
len=Length[B];
assoc=AssociationThread[B->Range[len]];
SymmetricBasis=SplitSymmetricBasis[B];
(*Seleccionar solo el sector deseado*)entries=SymmetricBasis[sector]/. {a_,b_,x_,y_}:>{{assoc[a],assoc[b]},x,y};
(*Construir base para el sector*)newBasis=Table[Module[{i1=e[[1,1]],i2=e[[1,2]],x=e[[2]],y=e[[3]]},ScaledCanonicalVectorij[i1,i2,x,y,len]],{e,entries}];
U1=SparseArray[newBasis];(*Matriz unitaria del sector*)Hsec=U1 . Horig . Transpose[U1];(*Hamiltoniano reducido*)Chop[Hsec,10^-10]  (*Eliminar componentes peque\[NTilde]as*)]




CreateHistogram[list_]:=Histogram[list,{Min[list],Max[list],(2*InterquartileRange[list])/(n^(1/3))},"PDF"]
KLevelSpacing[list_,order_]:=Module[{k=order,Sortedlist=Sort[list],n=Length[list]},
Table[(Sortedlist[[i+2*k]]-Sortedlist[[i+k]])/(Sortedlist[[i+k]]-Sortedlist[[i]]),{i,1,n-2*k}]]


End[]
EndPackage[]

