#
# InduceReduce.gi        The InduceReduce package         Jonathan Gruber
#
# Implementations
#
# Copyright (C) 2018     Jonathan Gruber
#
#############################################################################
##
## Initialize the options record with the default values.
## The meaning of the components is as follows:
##
## UseFlintLLL
##
## a boolean variable that tells the program to use experimental code to
## interface with flints implementation of the LLL algorithm.
## If false, GAPs builtin LLLReducedGramMat will be used.
##
## UsePcPresentation
##
## a boolean variable that tells the program to convert the p-parts of
## elementary subgroups to PcGroups before computing their irreducible
## characters. Otherwise their inherited subgroup
## representation will be used.
##
## PreComputePowMaps
##
## A boolean that indicates whether powermaps should be computed first.
## In this case the inverse mapping is also used to compute inner products.
##
## DoCyclicFirst
##
## a boolean variable that tells the program to induce the charaters of all
## cyclic subgroups before proceeding to the non-cyclic elementary subgroups
##
## DoCyclicLast
##
## a boolean variable that tells the program to induce the charaters of all
## cyclic subgroups before proceeding to the non-cyclic elementary subgroups
##
## LLLOffset
##
## an integer variable that tells the program to postpone the first LLL
## reduction until the characters of the LLLOffset-th elementary subgroup
## have been induced.
##
## DELTA
##
## a Float that specifies the parameter delta for the LLL reduction
##
## UseFiniteFields
##
## a boolean variable that tells the program to use finite field arithmetic
##
CTUngerDefaultOptions := rec(
	UseFlintLLL := true,
	UsePcPresentation := true,
	PreComputePowMaps := false,
	DoCyclicFirst := false,
	DoCyclicLast := false,
	LLLOffset := 0,
	Delta := 0.75,
	UseFiniteFields := false
);

#############################################################################
##
#V IndRed
##
## a record containing subfunctions of CharacterTableUnger
##
InstallValue( IndRed , rec(

#############################################################################
##
#F IndRed.Init( <G> )
##
## Initialization of the algorithm, returns a record which contains all
## relevant data concerning the group G
##
	Init:=function(G, Opt)
	local GR,i, TR, z, prime;
		TR:=IndRed.GroupTools(); # get group tools and reduce tools
	
		GR:=rec();
	
		GR.G:=G;   #  The group itself
		GR.C:=CharacterTable(G);  # the character table (without irreducible characters)
		GR.n:= Size(G);  # the order of G
		GR.classes:=ShallowCopy(ConjugacyClasses(GR.C)); # the conjugacy classes (mutable)
		GR.k:=Size(GR.classes); # number of conjugacy classes
		GR.classreps:=List(GR.classes,x->Representative(x)); # class representatives
		GR.orders:=List(GR.classreps, x->Order(x)); # element orders
			Info(InfoCTUnger, 1, "Induce/Restrict: group with ",
				Length(GR.orders), " conjugacy classes.");
			Info(InfoCTUnger, 2, "Induce/Restrict: orders of class reps: ",
				GR.orders);
		GR.perm:=Sortex(GR.orders,function(x,y) return y<x; end);
			# sort by decreasing order of representatives
		GR.classes:=Permuted(GR.classes,GR.perm); # adjust the order of classes and classreps
		GR.classreps:=Permuted(GR.classreps,GR.perm);
		GR.ccsizes:=List(GR.classes,x->Size(x)); # sizes of conjugacy classes
		GR.OrderPrimes:=List(GR.orders,x->Set(Factors(x)));
			# primes dividing the orders of class representatives
		GR.CentralizerPrimes:=List([1..GR.k],
			x-> Filtered( Set(Factors(GR.n/GR.ccsizes[x])) , y-> not y in GR.OrderPrimes[x] )  );
			# primes dividing the order of the centralizer, but not the order of the element
		GR.InduceCyclic:=ListWithIdenticalEntries(GR.k,true);
			# the characters of the corresponding cyclic groups still need to be induced
		GR.NumberOfCyclic:=GR.k; # number of cyclic groups whose characters need to be induced
		GR.NumberOfPrimes:=Sum(GR.CentralizerPrimes,x->Size(x));
			# number of primes for elementary subgroups not used so far
		GR.e:=Exponent(G); # exponent of G
		GR.ordersPos:=[];
			# create a non-dense list such that ordersPos[i] is the position of the first class of elements of order i
		for i in [1..GR.k] do
			if not IsBound(GR.ordersPos[GR.orders[i]]) then
				GR.ordersPos[GR.orders[i]]:=i;
			fi;
		od;
		GR.Ir:=[ Zero([1..GR.k])+1 ]; # initialize Irr(G) with the trivial character:
		GR.B:=[];  # B as an empty list
		GR.Gram:=[]; # Gram empty
		GR.IndexCyc:=0;
		# position of the cyclic group curretly used in the list of class representatives
		GR.m:=0; # positions from which on characters are not reduced so far
		GR.centralizers:=[];  # centralizers and powermaps computed so far
		GR.powermaps:=[];
		if Opt.PreComputePowMaps or Opt.DoCyclicFirst or Opt.UseFiniteFields then
			for i in [1..GR.k] do
				TR.PowMap(GR, i);
			od;
			GR.inverseClasses := List([1..GR.k], i -> GR.powermaps[i][GR.orders[i]]);
		fi;
		if Opt.UseFiniteFields then
			prime:=2*GR.n+1;

			while not IsPrimeInt(prime) do
				prime:=prime+GR.e;
			od;

			Info(InfoCTUnger, 1 ,"Choosing prime ", prime);

			GR.f := GF(prime);
			GR.omegaInt := PowerModInt(PrimitiveRootMod(prime),(prime-1)/GR.e,prime);
			GR.omega := GR.omegaInt * One(GR.f);
			GR.Ir:=[ Zero([1..GR.k])+One(GR.f) ];
			GR.ModularMap := DxGeneratePrimeCyclotomic(GR.e,GR.omega);
			GR.prime := prime;
			GR.normlimit := Int(prime/2);
			GR.ninv := One(GR.f)/GR.n;
		else
			GR.omega := E(GR.e);
			GR.ninv := 1/GR.n;
		fi;
		return GR;
	end ,

	#############################################################################
	##
	#F  DxLiftCharacter(<D>,<modChi>) . recalculate character in characteristic 0
	##
	LiftCharacter := function(GR, modular)
		local modularchi,chi,zeta,degree,sum,M,l,s,n,j,polynom,chipolynom,
			family,prime;

		Info(InfoCTUnger, 2, "Lifting modular character ", modular, ".\n");


		prime:=GR.prime;
		modularchi:=List(modular, Int);
		degree:=modularchi[GR.k];
		chi:=[];
		chi[GR.k] := degree;
		for j in [1..(GR.k)-1] do
			# FIXME
			# we need to compute the polynomial only for prime classes. Powers are
			# obtained by simply inserting powers in this polynomial
			family := GR.powermaps[j]{[2..GR.orders[j]]};
			l:=GR.orders[j];
			zeta:=E(l);
			polynom:=[degree,modularchi[j]];
			for n in [2..l-1] do
				s:=family[n];
				polynom[n+1]:=modularchi[s];
			od;
			chipolynom:=[];
			s:=0;
			sum:=degree;
			while sum>0 do
				M:=DxModularValuePol(polynom,
									PowerModInt(GR.omegaInt,-s*GR.e/l,prime),
									#PowerModInt(D.z,-s*D.irrexp/l,prime),
									prime)/l mod prime;
				Add(chipolynom,M);
				sum:=sum-M;
				s:=s+1;
			od;
			for n in [1..l-1] do
				s:=family[n];
				if not IsBound(chi[s]) then
					chi[s]:=ValuePol(chipolynom,zeta^n);
				fi;
			od;
		od;

		Info(InfoCTUnger, 2, "Lifted Character: ", chi, "\n");
		return chi;
	end,


#############################################################################
##
#F IndRed.GroupTools( )
##
## returns a record with different functions used for the algorithm
##
	GroupTools:=function()
	local TR;
		TR:=rec();

		## compute q-part of an integer n
		TR.pPart:= function(n,q)
		local res;
			res:=1;
			while RemInt(n,res)=0 do
				res:=res*q;
			od;
			return res/q;
		end;

		# find the position of the conjugacy class containing g in GR.classes ,
		# ord is the order of g
		TR.FindClass:= function(GR,h,ord)
		local j;
			for j in [GR.ordersPos[ord]..GR.k] do
				if not GR.orders[j]=ord then break; fi;	
				if h=GR.classreps[j] then 
				# first check if h actually equals one of the class representatives
					return j;
				fi;
			od;
			for j in [GR.ordersPos[ord]..GR.k] do
				if h in GR.classes[j] then
					return j;
				fi;
			od;
		end;

		## compute powermap of l-th class representative and add it to GR
		TR.PowMap:= function(GR,l)
		local h, res, i, ord;
			if GR.orders[l]=1 then GR.powermaps[l]:=[l]; fi; 
			res:=[TR.FindClass(GR,One(GR.G),1),l];
			h:=GR.classreps[l]^2;
			for i in [2..GR.orders[l]-1] do
				ord:=GR.orders[l]/GcdInt(GR.orders[l],i);
	
				Add(res , TR.FindClass(GR,h,ord));
			
				h:=h*GR.classreps[l];
			od;
			GR.powermaps[l]:=res;
		end;

		## Compute fusion of conjugacy classes of GR.Elementary to classes of G
		TR.ClassFusion:= function(GR)
		local i, j, res, ord, found;
			res:=[];
			for i in [1..GR.Elementary.k] do
				ord:=Order(GR.Elementary.classreps[i]);
				Add(res, TR.FindClass(GR,GR.Elementary.classreps[i],ord));
			od;
			return res;
		end;
	
		## Compute character table of cyclic group as matrix
		TR.Vandermonde:= function(GR)
		local i, j, M, omega;
			omega := GR.omega^(GR.e / GR.orders[GR.IndexCyc]);

			M:=NullMat(GR.orders[GR.IndexCyc],GR.orders[GR.IndexCyc]);
			for i in [0..GR.orders[GR.IndexCyc]-1] do
				for j in [0..GR.orders[GR.IndexCyc]-1] do
					M[i+1][j+1] := omega^(i*j);
				od;
			od;
			return M;
		end;

		# conjugacy class representatives of cyclic group corresponding to
		# the ordering of columns in Vandermonde
		TR.ClassesCyclic:=function(GR) 
		local res, h;
			res:=[Identity(GR.G)];
			h:=GR.Elementary.z;
			while not h=Identity(GR.G) do
				Add(res,h);
				h:=h*GR.Elementary.z;
			od;
			return res;
		end;

		return TR;

	end ,

#############################################################################
##
#F IndRed.ReduceTools( )
##
## returns a record with different functions required for lattice reduction
##
	# ReduceTools:=function()
	# local RedTR,ip, ipSparse;

	# 	RedTR:=rec();

	# 	# inner product of class functions x and y
		ip:= function(GR,x,y)
			local res;
			if IsBound(GR.inverseClasses) then
				res := Int(Sum([1..GR.k], i->x[i]*y[GR.inverseClasses[i]]*GR.ccsizes[i])*GR.ninv);
				if IsBound(GR.normlimit) then
					if res >= GR.normlimit then
						res := res - GR.prime;
					fi;
				fi;
				return res;
			else
				return Sum([1..GR.k], i->x[i]*ComplexConjugate(y[i])*GR.ccsizes[i])/GR.n;
			fi;
		end,
	
		# inner product of class functions x and y with few entries
		# one of them has all non-zero entries in GR.Elementary.FusedClasses
		ipSparse:= function(GR,x,y)
			local res;
			if IsBound(GR.inverseClasses) then
				res := Int(Sum(GR.Elementary.FusedClasses,
					i->x[i]*y[GR.inverseClasses[i]]*GR.ccsizes[i])*GR.ninv);
				if IsBound(GR.normlimit) then
					if res >= GR.normlimit then
						res := res - GR.prime;
					fi;
				fi;
				return res;
			else
				return Sum(GR.Elementary.FusedClasses,
					i->x[i]*ComplexConjugate(y[i])*GR.ccsizes[i])/GR.n;
			fi;
		end,
	
		# reduce new induced characters by the irreducibles found so far
		redsparse:=function(GR)
		local mat,pos;
			if GR.m+1>Size(GR.B) then return; fi;
			pos:=[GR.m+1..Size(GR.B)];
			mat:=List( pos , x->List(GR.Ir, y -> IndRed.ipSparse(GR,GR.B[x],y) ) );
			GR.B{pos}:=GR.B{pos}-mat*GR.Ir;
		end,
	
		# extend the gram matrix with the scalar products of the new characters
		GramMatrixGR:=function(GR)
		local i,j,mat,b;
			b:=Size(GR.B);
			mat:=NullMat(b,b);
			mat{[1..GR.m]}{[1..GR.m]}:=GR.Gram;
			for i in [GR.m+1..b] do
				for j in [1..i] do
					mat[i][j]:=IndRed.ip(GR,GR.B[i],GR.B[j]);
					mat[j][i]:=mat[i][j];
				od;
			od;
			GR.Gram:=mat;
		end,

	# 	return RedTR;
	# end ,

#############################################################################
##
#F IndRed.FindElementary( <GR> )
##
## finds the next elementary subgroup and puts it in the component Elementary
## of the record GR.
##
	FindElementary:=function(GR,Opt)
	local Elementary;
		Elementary:=rec();
		while true do
			GR.IndexCyc:=GR.IndexCyc+1; # run over the conjugacy classes
			if GR.IndexCyc>GR.k then GR.IndexCyc:=1; fi;
			if not IsEmpty(GR.CentralizerPrimes[GR.IndexCyc]) then
				# search for a prime that has not been used yet
				Elementary.isCyclic:=false;
				Elementary.z:=GR.classreps[GR.IndexCyc];
				Elementary.p:=Random(GR.CentralizerPrimes[GR.IndexCyc]);
					# choose a random prime for that conjugacy class
				RemoveSet(GR.CentralizerPrimes[GR.IndexCyc] , Elementary.p);
				GR.NumberOfPrimes:=GR.NumberOfPrimes-1;
				if not IsBound(GR.centralizers[GR.IndexCyc]) then
					# compute the centralizer, if not done already
					GR.centralizers[GR.IndexCyc]:=Centralizer(GR.G,Elementary.z);
				fi;
				Elementary.P:=SylowSubgroup(GR.centralizers[GR.IndexCyc],Elementary.p);
				GR.Elementary:=Elementary;
				break;
			elif GR.InduceCyclic[GR.IndexCyc] and
				(not Opt.DoCyclicLast or GR.NumberOfPrimes<=0) then
				# if there are no primes left for this class, induce from cyclic group
				# if Opt.DoCyclicLast=true: only induce from cyclic groups
				# when all primes have been used
				Elementary.isCyclic:=true;
				Elementary.z:=GR.classreps[GR.IndexCyc];
				GR.InduceCyclic[GR.IndexCyc]:=false;
				GR.NumberOfCyclic:=GR.NumberOfCyclic-1;
				GR.Elementary:=Elementary;
				break;
			fi;
		od;
		return;
	end ,

#############################################################################
##
#F IndRed.InitElementary( <GR>, <TR>, <Opt> )
##
## initialize data concerning the elementary subgroups and exclude its
## subgroups and their conjugates from the computations yet to come.
##
	InitElementary:=function(GR,TR,Opt)
		local i,j,i1,j1,p,temp,powermap;
		if GR.Elementary.isCyclic then
			GR.Elementary.n:=GR.orders[GR.IndexCyc]; # order of the elementary group
			GR.Elementary.k:=GR.Elementary.n; # number of conjugacy classes
			GR.Elementary.ccsizes:=ListWithIdenticalEntries(GR.Elementary.k,1); # class sizes
			if not IsBound(GR.powermaps[GR.IndexCyc]) then # if necessary, compute powermap
				TR.PowMap(GR,GR.IndexCyc);
			fi;
			powermap:=GR.powermaps[GR.IndexCyc];
			GR.Elementary.classfusion:=powermap; 
				# class fusion equals power map for cyclic group
			GR.Elementary.classreps:=TR.ClassesCyclic(GR); # class representatives
			GR.Elementary.XE:=TR.Vandermonde(GR); # character table
		else
			GR.Elementary.n:=GR.orders[GR.IndexCyc]*Size(GR.Elementary.P); # order
			if Opt.UsePcPresentation then
				GR.Elementary.pcIso := IsomorphismPcGroup(GR.Elementary.P);
				GR.Elementary.pcP := Image(GR.Elementary.pcIso);
				GR.Elementary.pcCtblP := CharacterTable(GR.Elementary.pcP);
				GR.Elementary.pcClassrepsP := List(ConjugacyClasses(GR.Elementary.pcCtblP), x->Representative(x));
				GR.Elementary.classrepsP := List(GR.Elementary.pcClassrepsP, x -> PreImagesRepresentative(GR.Elementary.pcIso, x));
				GR.Elementary.ccsizesP:=List(ConjugacyClasses(GR.Elementary.pcCtblP),x->Size(x));
				GR.Elementary.XP:=Irr(GR.Elementary.pcCtblP);
			else
				GR.Elementary.ctblP:=CharacterTable(GR.Elementary.P);
					# character table of p-group
				GR.Elementary.classrepsP:=List(ConjugacyClasses(GR.Elementary.ctblP),
					x->Representative(x));
					#classes of p-group in corresponding order
				GR.Elementary.ccsizesP:=List(ConjugacyClasses(GR.Elementary.ctblP),x->Size(x)); 				# class sizes p-group
				GR.Elementary.XP:=Irr(GR.Elementary.ctblP);
					# irreducible characters of the p-group
			fi;
			GR.Elementary.kP:=Size(GR.Elementary.classrepsP); # number of classes of p-group
			GR.Elementary.classrepsZ:=TR.ClassesCyclic(GR);
				# class representatives cyclic group
			GR.Elementary.kZ:=GR.orders[GR.IndexCyc]; # number of classes cyclic group
			if not IsBound(GR.powermaps[GR.IndexCyc]) then # if necessary, compute powermap
				TR.PowMap(GR,GR.IndexCyc);
			fi;
			powermap:=GR.powermaps[GR.IndexCyc];
			GR.Elementary.k:=GR.Elementary.kP*GR.Elementary.kZ;
				# number of classes elementary group
			GR.Elementary.classreps:=[]; # compute the class representatives
			GR.Elementary.ccsizes:=[]; # and the class sizes
			for i in [1..GR.Elementary.kP] do
				for j in [1..GR.Elementary.kZ] do
					Add(GR.Elementary.classreps ,
						GR.Elementary.classrepsP[i]*GR.Elementary.classrepsZ[j]);
					Add(GR.Elementary.ccsizes,GR.Elementary.ccsizesP[i]);
				od;
			od;
			GR.Elementary.classfusion:=TR.ClassFusion(GR);
				# compute the fusion of conjugacy classes of the elementary group
				# to conjugacy classes of G
			GR.Elementary.XZ:=TR.Vandermonde(GR); # character table of the cycic group
			if Opt.UseFiniteFields then
				GR.Elementary.XP := List(GR.Elementary.XP, chi -> List(chi, GR.ModularMap));
			fi;
			GR.Elementary.XE:=[]; # compute character table of the elementary group
			for i in [1..GR.Elementary.kP] do
				for j in [1..GR.Elementary.kZ] do
					temp:=[];
					for i1 in [1..GR.Elementary.kP] do
						for j1 in [1..GR.Elementary.kZ] do
							Add(temp,GR.Elementary.XP[i][i1]*GR.Elementary.XZ[j][j1]);
						od;
					od;
					Add(GR.Elementary.XE,temp);
				od;
			od;
		fi;
		GR.Elementary.FusedClasses:=Set(GR.Elementary.classfusion); 
			# positions of classes of G which contain classes of the elementary group
		for i in GR.Elementary.FusedClasses do # eliminate some elementary subgroups:
			if GR.InduceCyclic[i] then # the cyclic subgroups of the elementary groups
				GR.InduceCyclic[i]:=false;
				GR.NumberOfCyclic:=GR.NumberOfCyclic-1;
			fi;
			if GR.Elementary.isCyclic then
				# some elementary groups contained in the cyclic group
				for p in GR.CentralizerPrimes[i] do
					if TR.pPart(GR.n/GR.ccsizes[i],p)=TR.pPart(GR.orders[GR.IndexCyc],p) then
						RemoveSet(GR.CentralizerPrimes[i],p);
						GR.NumberOfPrimes:=GR.NumberOfPrimes-1;
					fi;
				od;
			fi;
		od;
		if not GR.Elementary.isCyclic then
			for i in Set(powermap) do 
				# elementary subgroups, where the cyclic part is generated
				# by a power of GR.Elementary.z and the p-groups coincide
				if TR.pPart(GR.n/GR.ccsizes[i],GR.Elementary.p) =
					TR.pPart(GR.n/GR.ccsizes[GR.IndexCyc],GR.Elementary.p) and
					GR.Elementary.p in GR.CentralizerPrimes[i] then
					RemoveSet(GR.CentralizerPrimes[i],GR.Elementary.p);
					GR.NumberOfPrimes:=GR.NumberOfPrimes-1;
				fi;
			od;
		fi;
		for i in [0..GR.orders[GR.IndexCyc]-1] do
			# derive the powermaps of powers of GR.Elementary.z
			if not IsBound(GR.powermaps[powermap[i+1]]) then
				GR.powermaps[powermap[i+1]]:=[];
				for j in [0..GR.orders[powermap[i+1]]-1] do
					Add(GR.powermaps[powermap[i+1]],
						powermap[ i*j mod GR.orders[GR.IndexCyc] +1 ]);
				od;
			fi;
		od;
		return;
	end ,

#############################################################################
##
#F IndRed.InduceCyc( <GR> )
##
## induce characters from all cyclic subgroups
##
	InduceCyc:=function(GR)
	local ords, inds;
		GR.Elementary:=rec();
		ords := ShallowCopy(OrdersClassRepresentatives(GR.C));
		inds := [1..NrConjugacyClasses(GR.C)];
		SortParallel(ords, inds, function(x,y) return x>y; end);
		GR.Elementary.FusedClasses := [1..GR.k];
		Append(GR.B, List(InducedCyclic(GR.C,inds,"all"), ch-> Permuted(ch, GR.perm)));
		return;
	end ,

#############################################################################
##
#F IndRed.Induce( <GR> )
##
## Induce all irreducible characters of the elementary subgroup in 
## GR.Elementary to the group GR.G and add them to GR.B
##
	Induce:=function(GR)
	local mat, i, j;
		mat:=NullMat(GR.Elementary.k,GR.k);
		for i in [1..GR.Elementary.k] do
			for j in [1..GR.Elementary.k] do
				mat[i][GR.Elementary.classfusion[j]]:=mat[i][GR.Elementary.classfusion[j]]+
					( (GR.n/GR.ccsizes[GR.Elementary.classfusion[j]]) / 
					(GR.Elementary.n/GR.Elementary.ccsizes[j])*GR.Elementary.XE[i][j] );
			od;
		od;
		mat:=Set(mat); # remove duplicates
		Append(GR.B,mat);
		return;
	end ,

#############################################################################
#
# IndRed.Reduce( <GR>, ><RedTR> )
#
# reduce the induced characters by the irreducibles found so far and do LLL
# lattice reduction
#
	Reduce:=function(GR,RedTR,Opt)
	local mat,temp,ind,I,i;
		IndRed.redsparse(GR); #reduce new characters by all irreducibles
		IndRed.GramMatrixGR(GR); # update the gram matrix
		if Opt.LLLOffset>0 then
			GR.m:=Size(GR.Gram);
		else
			if Opt.UseFlintLLL then
				temp:=LLLReducedGramMatFLINT(GR.Gram,Opt.Delta); # LLL reduction on gram matrix
			else
				temp:=LLLReducedGramMat(GR.Gram,Opt.Delta); # LLL reduction on gram matrix
			fi;
			GR.Gram:=temp.remainder;
			GR.B:=temp.transformation*GR.B;
			temp:=[];
			ind:=[];
			GR.m:=Size(GR.Gram);
			for i in [1..GR.m] do
				if Opt.UseFiniteFields then
					if GR.Gram[i][i] > GR.normlimit then
						Error("The chosen prime was too small. Tough luck.");
					fi;
				fi;
				if GR.Gram[i][i]=1 then
					Add(temp,i); # find positions of characters of norm 1
				else
					Add(ind,i);
				fi;
			od;
			if Size(temp)=0 then
				Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR.Ir),
					", dim = ", Length(GR.Ir)+Length(GR.Gram),
					", det(G) = ", DeterminantMat(GR.Gram));
				return;
			fi;
			I:=GR.B{temp}; # irreducible characters in B (up to sign)
			if Size(ind)=0 then
				Append(GR.Ir,I);
				GR.Gram:=[];
				Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR.Ir));
			return ;
			fi;
			for i in Reversed(temp) do
				Remove(GR.B,i); # remove irreducible characters from B
			od;
			mat:=GR.Gram{ind}{temp};
			GR.Gram:=GR.Gram{ind}{ind}-mat*TransposedMat(mat);
			GR.m:=Size(ind);
			GR.B:=GR.B-mat*I;
			Append(GR.Ir,Set(I)); # add the new irreducible characters to I
			Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR.Ir),
				", dim = ", Length(GR.Ir)+Length(GR.Gram),
				", det(G) = ", DeterminantMat(GR.Gram));
		fi;
		return;
	end ,

#############################################################################
##
#F IndRed.tinI( <GR> )
##
## adjusts the signs of the characters and permutes them to the right
## ordering of conjugacy classes
##
	tinI:=function(GR)
	local irr,i,perm;
		Print("Unit position: ", GR.ordersPos[1], "\n");

		if IsBound(GR.prime) then
			for i in [1..GR.k] do # adjust the signs
				if Int(GR.Ir[i][GR.ordersPos[1]]) > GR.normlimit then
					GR.Ir[i]:=-GR.Ir[i];
				fi;
			od;
			GR.Ir := List(GR.Ir, chi -> IndRed.LiftCharacter(GR, chi));
		else
			for i in [1..GR.k] do # adjust the signs
				if GR.Ir[i][GR.ordersPos[1]]<0 then
					GR.Ir[i]:=-GR.Ir[i];
				fi;
			od;
		fi;
		irr:=[];
		for i in [1..GR.k] do 
			# permute irreducible characters back to the order of classes in GR.C
			GR.Ir[i]:=Permuted(GR.Ir[i],Inverse(GR.perm));
			irr[i]:=Character(GR.C,GR.Ir[i]);
		od;
		GR.Ir:=irr;
		return;
	end ));
## end of the record IndRed #################################################

#############################################################################
##
#F InduceReduce(GR)
##
## computes the irreducible characters of the group GR.G (upto sign) and puts
## them in the component GR.Ir of GR, returns GR.Ir
##
InstallGlobalFunction( InduceReduce,
function(GR,Opt)
local TR, RedTR, H, ccsizesH, temp, it;

	TR:=IndRed.GroupTools(); # get group tools and reduce tools
	RedTR:= rec();


	if Opt.DoCyclicFirst = true then # if option Opt.DoCyclicFirst is set,
		# induce from all cyclic groups first and reduce

		Info(InfoCTUnger, 2, "Induce: from cyclic subgroups");

		IndRed.InduceCyc(GR);
		IndRed.Reduce(GR,RedTR,Opt);
		GR.InduceCyclic:=ListWithIdenticalEntries(GR.k,false);
		GR.NumberOfCyclic:=0;
	fi;

	if Size(GR.Ir)>=GR.k then 
		return GR.Ir;
	fi;


	while Size(GR.Ir)<GR.k and not (GR.NumberOfPrimes=0 and GR.NumberOfCyclic=0) do
		# if number of irr. characters=number of conjugacy classes,
		# all irr. characters have been found.
	
		Opt.LLLOffset:=Opt.LLLOffset-1;
		# Opt.LLLOffset postpones the first LLL lattice reduction
	
		IndRed.FindElementary(GR,Opt); # find elementary subgroup
	
		IndRed.InitElementary(GR,TR,Opt); # determine information needed about elementary subgroup
	
	
		Info(InfoCTUnger, 1, "Induce/Restrict: Trying [|Z|, |P|, k(E)] = ",
			[ GR.orders[GR.IndexCyc], 
				GR.Elementary.n/GR.orders[GR.IndexCyc], GR.Elementary.k ]);

		IndRed.Induce(GR); # append induced characters to GR.B
	
		IndRed.Reduce(GR,RedTR,Opt); # reduce GR.B by GR.Ir and do lattice reduction
	od;

	if Size(GR.Ir)>=GR.k then
		return GR.Ir;
	else
		Opt.Delta:=1; 
		# if there are still characters missing, do on last LLL reduction with Opt.Delta:=1
		Opt.LLLOffset:=0;
		IndRed.Reduce(GR,RedTR,Opt);
		return GR.Ir;
	fi;

end );

#############################################################################
##
#F CharacterTableUnger( <G> )
## 
## Computes the character table of a finite group using Unger's algorithm
##
InstallGlobalFunction( CharacterTableUnger,
function(G, Options...)
local GR, Opt, T;

	Opt:=ShallowCopy(CTUngerDefaultOptions);
	if Length(Options)>0 and IsRecord(Options[1]) then
		if IsBound(Options[1].DoCyclicFirst) and IsBool(Options[1].DoCyclicFirst) then
			Opt.DoCyclicFirst:=Options[1].DoCyclicFirst;
		fi;
		if IsBound(Options[1].DoCyclicLast) and IsBool(Options[1].DoCyclicLast) then
			Opt.DoCyclicLast:=Options[1].DoCyclicLast;
		fi;
		if IsBound(Options[1].LLLOffset) and IsInt(Options[1].LLLOffset) then
			Opt.LLLOffset:=Options[1].LLLOffset;
		fi;
		if IsBound(Options[1].Delta) and IsRat(Options[1].Delta) 
			and 1/4<Options[1].Delta and Options[1].Delta<=1 then
			Opt.Delta:=Options[1].Delta;
		fi;
	fi;

	GR:=IndRed.Init(G, Opt);

	InduceReduce(GR,Opt); # do the induce-reduce algorithm

	IndRed.tinI(GR);
	
	# convert the result into a chacacter table
	T:=rec();
	T.UnderlyingCharacteristic:=0;
	T.NrConjugacyClasses:=GR.k;
	T.ConjugacyClasses:=ConjugacyClasses(GR.C);
	T.Size:=GR.n;
	T.OrdersClassRepresentatives:=Permuted(GR.orders,Inverse(GR.perm));
	T.SizesCentralizers:=Permuted(List(GR.ccsizes,x->GR.n/x),Inverse(GR.perm));
	T.Irr:=GR.Ir;
	T.UnderlyingGroup:=UnderlyingGroup(GR.C);
	T.IdentificationOfConjugacyClasses:=IdentificationOfConjugacyClasses(GR.C);
	T.InfoText:="Computed using Unger's induce-reduce algorithm";
	ConvertToCharacterTableNC(T);
	return T;
end );
