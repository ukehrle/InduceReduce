InstallGlobalFunction( InduceReduceParallel,
function(GR)
	local H, ccsizesH, temp, it, Elementary, tmp, p, i, Opt, wf, cnt, chars, fusedClasses, job;

	Opt := IRGetOptions();

    if Opt.IRDoCyclicFirst = true then
        Info(InfoCTUnger, 1, "Induce: from cyclic subgroups");
        IRDoCyclics(GR);
    fi;

	if Size(GR!.Ir)>=GR!.k then
		return GR!.Ir;
	fi;

    wf := List([1..Opt.IRParallelGroups], i -> BackgroundJobByFork(function(el) return IRHandleElementary(GR, el); end, fail ));

	while true  do
        chars := [];
        fusedClasses := [];

        cnt := 0;

        for job in wf do
            if not IsIdle(job) = true then
              continue;
            fi;

            tmp := Pickup(job);
            if tmp <> false then
                Info(InfoCTUnger, 2, "IR/Parallel: Added new characters");
                UniteSet(chars, tmp.characters);
                UniteSet(fusedClasses, tmp.fusedClasses);
                cnt := cnt + 1;
            fi;

            Elementary := IRFindGroup(GR); # find elementary subgroup
            if Elementary <> fail then
                Info(InfoCTUnger, 2, "IR/Parallel: Submitting new job.");
                Submit(job, [Elementary]);
            else
                # no more groups, kill worker.
                Kill(job);
            fi;
            break;
        od;

        if cnt = 0 then
          continue;
        fi;

		Opt.IRLLLOffset := Opt.IRLLLOffset - cnt;

        GR!.B := Concatenation(GR!.B, chars);

        Info(InfoCTUnger, 2, "IR/Parallel: entering reduce stage.");
        Info(InfoCTUnger, 2, GR);
		IRReduce(GR, fusedClasses : Opt); # reduce GR!.B by GR!.Ir and do lattice reduction
        Info(InfoCTUnger, 2, "IR/Parallel: exiting reduce stage.");

        if IRIsComplete(GR) then
            break;
        fi;
	od;

    for job in wf do
      Kill(job);
    od;

	if Size(GR!.Ir)>=GR!.k then
		return GR!.Ir;
	else
		# if there are still characters missing, try one last LLL reduction with Opt.Delta = 1
		Opt.IRDelta:=1;
		Opt.IRLLLOffset := 0;
		IRReduce(GR, tmp.fusedClasses : Opt);
		if Size(GR!.Ir)>=GR!.k then
			return GR!.Ir;
		fi;

		OnBreak := function() end;
		OnBreakMessage := function() end;
		Error("**********************************************************",
			  "\n\n",
			  "Run out of subgroups and still did not find all irreducible ",
			  "characters!\n",
			  "They are guaranteed to be integer linear combinations of\n",
			  "GR!.Ir and GR!.B, but lattice reduction did not find them.\n",
			  "You now have the opportunity to alter GR!.Ir in the REPL.\n",
			  "Type \"return;\" when you're done.\n\n");

		return GR!.Ir;
	fi;

end );
