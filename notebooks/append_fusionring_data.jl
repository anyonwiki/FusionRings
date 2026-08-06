using Revise, Oscar, Accessors, JSON, FusionRings, ProgressMeter

datadir = "/home/gert/Projects/FusionRings.jl/src/data/"

unique_rings =
  try
    if length(unique_rings) == 25138
      unique_rings
    else
      import_rings("/home/gert/Projects/FusionRings.jl/src/data/unique_frl.json")
    end
  catch e
    import_rings("/home/gert/Projects/FusionRings.jl/src/data/unique_frl.json")
  end


#updated_rings = FusionRing[]

#@showprogress dt=1 desc="Updating ring info" for fr in unique_rings
#  push!(updated_rings,auto_complete_missing_info(fr))
#end
#

#updated_rings2 = updated_rings
#@showprogress for i in 1:length(updated_rings)
#  subrings = sub_fusion_rings(updated_rings[i])
#  if isempty(subrings)
#    nr = updated_rings2[i]
#    nr = @set nr.sub_fusion_rings = Tuple{Vector{Int64},String}[]
#    updated_rings2[i] = nr
#  else
#    sfrstr = [ ( t[1], uuid(t[2]) ) for t in subrings ]
#    nr = updated_rings2[i]
#    nr = @set nr.sub_fusion_rings = sfrstr
#    updated_rings2[i] = nr
#  end
#end


#export_rings(joinpath(datadir,"updated_unique_frl.json"),updated_rings2)

#tabs = Oscar.load("src/data/unique_rk10_multtabs.mrdi")


#r10rings = FusionRing[]

function fix_central_series(fr)
  ucs = upper_central_series(fr;force_compute=true)

  ucss = [ (t[1],uuid(t[2])) for t in ucs ]

  change_properties(
    fr,
    :upper_central_series => ucss
  )
end


function fix_sfr(fr)
  sfr = sub_fusion_rings(fr;force_compute=true)

  sfrs = [ (t[1],uuid(t[2])) for t in sfr ]

  change_properties(
    fr,
    :sub_fusion_rings => sfrs
  )
end

function fix_r(fr)
  r = decompositions(fr;force_compute=true)

  sr = Dict{String,Any}( "tensor_product" => [ [ uuid(ring) for ring in decomp ] for decomp in r ] )

  change_properties(
    fr,
    :realizations => sr
  )
end
#@showprogress dt=1 desc="Updating ring info" for mt in tabs
#  push!(r10rings,fusion_ring(mt))
#end


function fr_fn(fr::FusionRing)
  c  = string.(anyonwiki_code(fr))
  us = fill("_", 3)
  return stringriffle(c, us) * ".json"
end

function export_new_rings(v::Vector{Int64})
  return export_rings(joinpath(datadir, "newrings.json"), import_and_fix.(frl[v]))
end
