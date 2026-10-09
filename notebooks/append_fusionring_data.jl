using Revise, Oscar, Accessors, JSON, FusionRings, ProgressMeter, ProfileView

# datadir = "/home/gert/Projects/FusionRings.jl/src/data/"

function update_rings( dir::String, frs::Vector{FusionRing} )
  bm  = group_by( multiplicity, frs )
  kys = collect(keys(bm))

  @showprogress for k in kys
    fn = dir*"proposals/AlgebraicStructures/FusionRings/"*"mult$k"*".json"
    export_rings( fn, bm[k] )
  end
end

function is_mrdi_file_name( fn::String )::Bool
  fn[end-4:end] == ".mrdi"
end

# qqb_fns = filter( is_mrdi_file_name, readdir(datadir*"split_number_data") );

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

texrepinfo = """
    JSON dictionary of LaTeX representations of algebraic numbers.

    The following conventions are used in explaining the values of the dictionary:
    * A qqb_id is a string that uniquely describes an algebraic number. It is formatted as a list of n integers a0, ..., an, separated by underscores, followed by a double underscore, followed by an integer i: "a0_a1_..._an__i". Here a0 to an are the coefficients of a polynomial a0 + a1*x + ... + an*x^n and i denotes the i'th root of that polynomial. The indexing of roots of the polynomial takes the real roots first, in increasing order. Then come complex conjugate pairs of roots, sorted first by increasing real part and second by increasing complex part.

    The "data" field maps qqb_id strings to dictionaries containing LaTeX representations of the algebraic number represented by the qqb_id or an empty string if we have no such representation available. Note that an empty string does not mean there exists no such representation. There are 3 representations available at the moment:
      * "radical": a representation of the number using a combination of sums products of rational powers of integers
      * "cyclo": a representation of the number using cyclotomic numbers where \\zeta_{n} stands for exp(2πi/n)
      * "general": a representation of the number as "[minpoly]_n" where minpoly is the minimal polynomial and n is the index of the root using the same indexing as for the root number of its qqb_id
"""

function update_algebraic_numbers(dir::String)
  vals = Oscar.load( dir*"SupportingData/AlgebraicNumbers/qqb_vals.mrdi");
  ids  = Oscar.load( dir*"SupportingData/AlgebraicNumbers/qqb_ids.mrdi");

  length(vals) != length(ids) && error("length(vals) != length(ids)")

  # create partitions of tuples of ids and vals of
  # size 10000
  ivp = collect(Iterators.partition( zip(ids,vals), 10000 ))

  @showprogress for i in 1:length(ivp)
    d = Dict{String,Any}(
      "info"  => texrepinfo,
      "order" => first.(ivp[i]),
      "data"  => Dict( t[1] => FusionRings.tex_reps(t[2],t[1]) for t in ivp[i]  )
    )
    FusionRings.write_json(dir*"proposals/SupportingData/AlgebraicNumbers/qqb_tex_reps_$i.json", d)
  end
  
end

function update_ring_info( rings )
  newrings = FusionRing[]

    sftw = Dict{String,Any}(
      "all" =>
          Dict(
              "name"      => "FusionRings",
              "version"   => "0.2.5",
              "swhid"     => "swh:1:rev:0dc2a749f029b30ca0e4c132d63432b4031ca68b;origin=https://github.com/anyonwiki/FusionRings;visit=swh:1:snp:efd4aefa347ea7a7b91d0de2da9a41ad47ffea09",
              "doi"       => nothing,
              "url"       => "https://github.com/anyonwiki/FusionRings"
          )
      )
  
  @showprogress for r in frl 
    if !is_commutative(r)
      push!(
        newrings,
        change_properties(
          sort(r),
          Dict(
            :software=>sftw,
            :formal_codegrees=>missing
          )
        )
      )
    else
      push!(
        newrings,
        change_properties(
          sort(r),
          :software=>sftw
        )
      )
    end
  end

  newrings
end


function correct_characters(ring)
  !is_commutative(ring) && return true

  ch = ComplexF64.(Matrix(characters(ring)))

  mt = multiplication_table(ring)

  mats = [ mt[i,:,:] for i in 1:rank(ring) ]
  
  action(mat) = ch*ComplexF64.(mat)*ch'

  sum(norm.(action.(mats)))
  
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
