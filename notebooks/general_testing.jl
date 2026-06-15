using Revise, Pkg, Oscar
Pkg.activate("/home/gert/Projects/FusionRings.jl/")
using FusionRings

function multtabcode(mt::Array{Int, 3}, mult::Int)::ZZRingElem
  m = mult + 1
  digits_base_m = [mt[i] for i in eachindex(mt)]
  digits_string = join(string.(digits_base_m))

  return ZZ(parse(BigInt, digits_string; base = m))
end

multtabcode(multiplication_table(frl[30]), 1)
