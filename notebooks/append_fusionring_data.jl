using Pkg;
Pkg.activate("/home/gert/Projects/FusionRings.jl")
using Revise, Oscar, Accessors, JSON, FusionRings, ProgressMeter

function to_qqb_arr(mat)
  return [mat[i] for i in eachindex(mat)]
end

# TODO: make sure numeric characters match the symbolic ones

function import_and_fix( )
  fn = "/home/gert/Projects/FusionRings.jl/src/data/fusionrings.json"

  js = JSON.parsefile(fn)

  # Mult tables were incorrect
  mt = multiplication_table(fr)



  return fusion_ring(
    mt;
    names                               = nfromjs(js),
    texnames                            = tnames,
    barcode                             = bcfromjs(js),
    anyonwiki_code                      = fcfromjs(js),
    characters                          = chfromjs(js),
    sub_fusion_rings                    = sfrfromjs(js),
    projective_SL2Z_reps                = psrfromjs(js),
    frobenius_perron_dimension          = fpdfromjs(js),
    frobenius_perron_dimensions         = fpdsfromjs(js),
    tensor_product_decompositions       = tpdfromjs(js),
    numeric_characters                  = nchfromjs(js),
    numeric_projective_SL2Z_reps        = npsrfromjs(js),
    numeric_frobenius_perron_dimension  = nfpdfromjs(js),
    numeric_frobenius_perron_dimensions = nfpdsfromjs(js),
    has_categories_with_props           = catswithprops,
    categorifications                   = ctsfromjs(js),
    references                          = js["references"],
    software                            = js["software"],
    comments                            = js["comments"],
  )
end

function fix_characters(fr::FusionRing)
  ch = fr.characters
  ch === missing && return fr

  @set fr.characters = string.(permutedims(ch))
end

function fr_fn(fr::FusionRing)
  c  = string.(anyonwiki_code(fr))
  us = fill("_", 3)
  return stringriffle(c, us) * ".json"
end

function export_new_rings(v::Vector{Int64})
  datadir = "/home/gert/Projects/FusionRings.jl/src/data/FusionRings/"
  return export_rings(joinpath(datadir, "newrings.json"), import_and_fix.(frl[v]))
end

macro timeout(seconds, expr_to_run, expr_when_fails = nothing)
  quote
    timer = Channel{Timer}(1)
    tsk = @task begin
      x = $(esc(expr_to_run))
      close(take!(timer))
      x
    end
    schedule(tsk)

    put!(
      timer,
      Timer($(esc(seconds))) do timer
        return istaskdone(tsk) || schedule(tsk, InterruptException(); error = true)
      end,
    )

    try
      fetch(tsk)
    catch _
      $(esc(expr_when_fails))
    end
  end
end
