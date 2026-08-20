using Revise, Oscar, Accessors, JSON, FusionRings, ProgressMeter, ProfileView

datadir = "/home/gert/Projects/FusionRings.jl/src/data/"



function update_rings( frs::Vector{FusionRing} )
  bm  = group_by( multiplicity, frs )
  kys = collect(keys(bm))

  @showprogress for k in kys
    fn = datadir*"FusionRings/"*"mult$k"*".json"
    export_rings( fn, bm[k] )
  end
end

function is_mrdi_file_name( fn::String )::Bool
  fn[end-4:end] == ".mrdi"
end

qqb_fns = filter( is_mrdi_file_name, readdir(datadir*"split_number_data") );

function qqb_id_from_fn(fn::String)::String
  fn[1:end-5]
end

function update_algebraic_numbers()

  fns = qqb_fns

  qqb_ids  = @showprogress desc = "Loading qqb ids"  [ qqb_id_from_fn(fn) for fn in fns ]
  Oscar.save(datadir * "AlgebraicNumbers/qqb_ids.mrdi",qqb_ids)

  function fn_to_path(fn::String)::String
    datadir*"split_number_data/"*fn
  end

  qqb_vals = @showprogress desc = "Loading qqb nums" [ Oscar.load(fn_to_path(fn)) for fn in fns ]
  Oscar.save(datadir * "AlgebraicNumbers/qqb_vals.mrdi",qqb_vals)

  function idandrep(i::Int64)::Tuple{String,Dict{String,String}}
    id = qqb_id_from_fn(fns[i])
    ( id, FusionRings.tex_reps(qqb_vals[i],id) )
  end

  qqb_tex_reps_arr = @showprogress desc = "Computing tex reps" [ idandrep(i) for i in 1:length(fns) ]
  trd = Dict( t[1] => t[2] for t in qqb_tex_reps_arr )

  FusionRings.write_json(datadir * "AlgebraicNumbers/qqb_tex_reps.mrdi", trd )
end

# update the categorifiability info
#catdata = JSON.parsefile("/home/gert/Projects/Lyctr.jl/src/data/MultFreeFusionCategories.json")["data"]
#=
# group the cats by fusionring
cats = collect(values(catdata))
groupedcats = group_by(c -> c["fusion_ring"] ,cats)

newfrl = frl;

uptorank7ids = uuid.( filter( r -> multiplicity(r)==1 && rank(r)<8, frl ) )

function update_ring(ring::FusionRing)
  id = uuid(ring)
  if haskey(groupedcats,id)
    cats = groupedcats[id]

    props = props_from_cats(cats)

    change_properties(
      ring,
      Dict(
      :has_categories_with_props => props,
      :categorifications => [ c["uuid"] for c in cats ],
      :software => Dict(
        "all_gradings" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "formal_codegrees" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "all_other_data" => ["https://doi.org/10.5281/zenodo.10686859"]
      ),
      :references => Dict(
        "all" => ["https://doi.org/10.1063/5.0148848"]
      )
      )
    )
  elseif id ∈ uptorank7ids
    change_properties(ring,
    Dict(
      :has_categories_with_props =>
        Dict(
          "Fusion"    => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
          "Pivotal"   => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
          "Unitary"   => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
          "Spherical" => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
          "Braided"   => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
          "Modular"   => [ false, "Computer:symbolic", "Anyonica v0.9.2 found no data" ],
        ),
      :categorifications => String[],
      :software => Dict(
        "all_gradings" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "formal_codegrees" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "all_other_data" => ["https://doi.org/10.5281/zenodo.10686859"]
      ),
      :references => Dict(
        "all" => ["https://doi.org/10.1063/5.0148848"]
      )
    )
    )
  else
    change_properties(ring,
    Dict(
      :software => Dict(
        "all_gradings" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "formal_codegrees" => ["https://github.com/anyonwiki/FusionRings@0.2.1"],
        "all_other_data" => ["https://doi.org/10.5281/zenodo.10686859"]
      ),
      :references => Dict(
        "all" => ["https://doi.org/10.1063/5.0148848"]
      )
    )
    )
  end
end

function props_from_cats(cats)
  Dict(
    "Fusion"    => [ true, "Computer:symbolic", "See categorifications" ],
    "Pivotal"   => [ true, "Computer:symbolic", "See categorifications" ],
    "Unitary"   => [ any([ c["is_unitary"] for c in cats ] ), "Computer:symbolic", "See categorifications" ],
    "Spherical" => [ any([ c["is_spherical"] for c in cats ] ), "Computer:symbolic", "See categorifications" ],
    "Braided"   => [ any([ c["is_braided"] for c in cats ] ), "Computer:symbolic", "See categorifications" ],
    "Modular"   => [ any([ c["is_modular"] for c in cats ] ), "Computer:symbolic", "See categorifications" ]
  )
end
#
#
=#

# adding names to new rings

# psu2k

# su2k

# grouprings

## zn rings

# rep fusion rings

# combinations & extensions of rings

# tensor products of rings

# HI

# TY

# adjoint rings
