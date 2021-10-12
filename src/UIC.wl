(* ::Package:: *)

BeginPackage["UIC`"];


GeneralizedLambertian::useage = "GeneralizedLambertian[I0,r,\[Theta],s] is a Generalized Lambertian model with inverse square law(\!\(\*FractionBox[\(I0\), SuperscriptBox[\(r\), \(2\)]]\)\!\(\*SuperscriptBox[\(cos\), \(s\)]\)(\[Theta])) for given I0: irradiation flux at r=1, r: distance bewtween source and target, s: Lambertian number, \[Theta]: angle between source direction and target. "


GeneralizedLamberXY::useage ="GeneralizedLambertian[I0,r,\[Theta],s] is a Generalized Lambertian model with vertical 2D plane, x: source location below plane, t: target location above plane, h: distance between two plane."


GeneralizedLamberXYPhase::useage  "GeneralizedLamberXYPhase[I0,x,t,h,s,\[Phi]] is a phase effect added generalized Lambertian 2D model. \[Phi]>0 outward rotation of source from center to boundary."


LamberIntensityCenter::useage ="LamberIntensityCenter[I0,x,h,s]: Center irradiance intensity of two s-generalized Lambertian light sources at x, -x location for h distance target plane"


LamberIntensityBoundary::useage ="Boundary irradiance intensity of two s-generalized Lambertian light sources at x, -x location for h distance target plane"


LamberD::useage ="LamberD[a,d,s]: Normalized boundary, center irradiance difference of two s-generalized Lambertian light sources at d, -d location"


LamberFindCorrespondingPointQR::useage ="LamberFindCorrespondingPointQR[dp,de,a,h,s]: Find corresponding points in Q-R region for point in such region using normalized LamberD function"


LamberFindde::useage ="LamberFindde[a,h,s]: using normalized LamberD function"


LamberFinddm::useage ="LamberFinddm[de,a,h,s]: using LamberFindCorrespondingPointQR function"


LamberFinddmApprox::useage ="LamberFinddmApprox[s]:"


LamberFinddeApprox::useage ="LamberFinddeApprox[s,a]:"


LamberDi::useage ="LamberDi[I0, x, h, w,s]:"


LamberFindCorrespondingxPointQR::useage ="LamberFindCorrespondingxPointQR[x,xe,a,h,s]:"


LamberFindxe::useage ="LamberFindxe[w_,h_,s_]:"


LamberFindxm::useage ="LamberFindxm[xe_,w_,h_,s_]:"


LamberFindxeApprox::useage ="LamberFindxmApprox[s_,h_]:"


LamberFindxmApprox::useage ="LamberFindxeApprox[s_,h_,w_]:"


xeApproxEsti::useage="xeApproxEsti[s,a]"


xmApproxEsti::useage="xmApproxEsti[s,a]"


ESCCoefficientLinear::useage ="ESCCoefficientLinear[n_,s_]:"


ESCCoefficientRectangular::useage ="ESCCoefficientRectangular[n_,m_,s_]:"


ESCArrayLinearN::useage ="ESCArrayLinearN[n_,s_,h_,half_]:"


ESCArrayRectangularNM::useage ="ESCArrayRectangularNM[n_,m_,s_,h_,half_]:"


(*ArrayIntensity::useage ="ArrayIntensity[x_,t_,s_,h_,w_]:"*)


Begin["`Private`"]


GeneralizedLambertian[I0_,r_,\[Theta]_,s_]:= I0/r^2*Power[Cos[\[Theta]],s]


GeneralizedLamberXY[I0_,x_,t_,h_,s_]:=I0*h^s/Power[(h^2+(t-x)^2),s/2+1]


GeneralizedLamberXYPhase[I0_,x_,t_,h_,s_,\[Phi]_]:= I0* 1/Power[(h^2+(t-x)^2),s/2+1]*(h*Cos[\[Phi]]+(t-x)*Sin[\[Phi]])^s


LamberIntensityCenter[I0_,x_,h_,s_]:= GeneralizedLamberXY[I0,x,0,h,s] + GeneralizedLamberXY[I0,-x,0,h,s] 


LamberIntensityBoundary[I0_,x_,h_,w_,s_]:=GeneralizedLamberXY[I0,x,w/2,h,s] + GeneralizedLamberXY[I0,-x,w/2,h,s] 


LamberD[a_,d_,s_]:= Power[(1+((1/2)*a + d)^2),-((s+2)/2)] + Power[(1+((1/2)*a - d)^2),-((s+2)/2)]-2*Power[(1+(d)^2),-((s+2)/2)]


LamberFindCorrespondingPointQR[dp_,de_,a_,h_,s_]:=
Module[{dc, l},
		l = LamberD[a,dp,s];
		dc = t/. FindRoot[LamberD[a,t,s]-l,{t, de/2+a/4}]
	Return[dc];
];


LamberFindde[a_,s_]:= d/. FindRoot[LamberD[a,d,s],{d,a/4}]


LamberFinddm[de_,a_,s_]:= 
Module[{d},
    d = x/. FindRoot[LamberD[a,x,s]+LamberD[a,a/2,s],{x,de*0.5}];
    Return[d];
];


LamberFinddmApprox[s_]:= Sqrt[Power[2,2/(s+2)]-1]


LamberFinddeApprox[s_,a_]:=1/6*LamberFinddmApprox[s] + 1/4*a


LamberDi[I0_, x_, h_, w_,s_] := (I0/h^2)*LamberD[w/h,x/h,s]


LamberFindCorrespondingxPointQR[x_,xe_,a_,h_,s_]:=h*LamberFindCorrespondingPointQR[x/h,xe/h,a,h,s];


LamberFindxe[w_,h_,s_]:= h* LamberFindde[w/h,s];


LamberFindxm[xe_,w_,h_,s_]:= h*LamberFinddm[xe/h,w/h,s];


LamberFindxmApprox[s_,h_]:=h*LamberFinddmApprox[s];


LamberFindxeApprox[s_,h_,w_]:=h*LamberFinddeApprox[s,w/h];


xeApproxEsti[s_,a_]:=Module[{x,l}, l=LamberFindde[a,s]; x = Abs[N[LamberFinddeApprox[s,a]]-l]/l *100; Return[x]];


xmApproxEsti[s_,a_]:=Module[{x,l}, de = LamberFindde[a,s]; l=LamberFinddm[de,a,s]; x = Abs[N[LamberFinddmApprox[s]]-l]/l *100; Return[x]];


ESCCoefficientLinear[n_,s_]:=
Module[{d,f,x0=0.6}, 
	If[Mod[s,2]<1,x0=0.8,x=0.6];
	f[x_]:=Sum[Power[(n+1-2*i)^2*(x/2)^2+1,-((s+6)/2)]*(1-(s+3)*(n+1-2*i)^2 *(x/2)^2),{i,n}];
	If[EvenQ[n],
		If[n==2,d=Power[4/(s+3),1/2], d =x/. FindRoot[f[x],{x,x0}];];,
		d =x/. FindMinimum[{f[x], 0<= x<= 1},{x,x0}][[2]];
	];
	Return[ d]
];



ESCCoefficientRectangular[n_,m_,s_]:=
Module[{d},
	f[x_]:= Sum[Power[((n+1-2*i)^2+(m+1-2*j)^2)*(x^2/4)+1, -(s+6)/2] *(1-((s+3)*(n+1-2*i)^2-(m+1-2*j)^2)*(x^2/4)),{j,m},{i,n}];
	if[n==m && n==2, 
		Return[Sqrt[4/(s+2)]]; v
		,
		(*Root try -> Minimize*)
         Quiet[
             Check[
                d = x/.FindRoot[f[x],{x,0.5}];
                , 
                d= x/.FindMinimum[{f[x], 0<= x<= 1},{x,0.5}][[2]];
            ];
        ];
	];
	Return[d];
];



ESCArrayLinearN[n_,s_,h_,half_]:=
Module[{},
	d= h* ESCCoefficientLinear[n,s];
	L = {n,d, d*(n-1)};
	If[half ==0,
		AppendTo[L, Table[(-(n-1)/2+i-1)*d,{i,n}]];
		,
		If[ EvenQ[n],
			AppendTo[L, Table[(1/2+i)*d,{i,n/2}]];
			,
			AppendTo[L, Table[(i-1)*d,{i,(n+1)/2}]];
			];
		];
		Return[L];
];



ESCArrayRectangularNM[n_,m_,s_,h_]:=
Module[{},
	d = h* ESCCoefficientRectangular[n,m,s];
	wx = (m-1)*d;
	wy = (n-1)*d;
	xarray = Table[];
	yarray = Table[];
	Return[Table[xarray,yarray]]
];


(*
ArrayIntensity[arr_,t_,s_,h_,w_]:=
Module[{n,i},
	m=Length[arr[1]];
	n=Length[arr[2]];
	i = N[];
	If[];
	Return[i];
];
*)


End[];


EndPackage[];
