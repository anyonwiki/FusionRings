# ┌────────────────────────────────────────────────────────────────────────────┐
# │                              Pretty printing                               │
# └────────────────────────────────────────────────────────────────────────────┘

function transform_integer(i::Int, dict::Dict)
  dgts = digits(i; base = 10)
  return join(reverse([dict[j] for j in dgts]))
end

bold_integer(i::Int)::String = transform_integer(i, bold_digits_dict)

bold_digits_dict = Dict(
  0 => "𝟎",
  1 => "𝟏",
  2 => "𝟐",
  3 => "𝟑",
  4 => "𝟒",
  5 => "𝟓",
  6 => "𝟔",
  7 => "𝟕",
  8 => "𝟖",
  9 => "𝟗",
)

subscript_integer(i::Int)::String = transform_integer(i, subs_digits_dict)

subs_digits_dict = Dict(
  0 => "₀",
  1 => "₁",
  2 => "₂",
  3 => "₃",
  4 => "₄",
  5 => "₅",
  6 => "₆",
  7 => "₇",
  8 => "₈",
  9 => "₉",
)

function superscript_integer(i::Int)
  return transform_integer(i, sup_digits_dict)
end

sup_digits_dict = Dict(
  0 => "⁰",
  1 => "¹",
  2 => "²",
  3 => "³",
  4 => "⁴",
  5 => "⁵",
  6 => "⁶",
  7 => "⁷",
  8 => "⁸",
  9 => "⁹",
)

# Formatting of fusion rings
# Names are grouped per category. The user can change the order of importance preference of which
# category is used to display the names of all rings

global _naming_priority_order = [ "quantum_group_like", "group_like", "physics", "miscellaneous" ]
global _naming_allow_exceptions
global _naming_exceptions = Dict(
  [ 2, 1, 0, 2 ] => "Fibonacci",
  [ 3, 1, 0, 1 ] => "Ising"
)

export set_naming_priority!

function set_naming_priority!(v::Vector{String})
  v ⊈ _naming_priority_order && error("The priorities given for naming should be a subset of [ \"quantum_group\", \"group_like\", \"physics\", \"miscellaneous\" ] ")
  np = v
  for s in _naming_priority_order
    s ∉ np && push!(np,s)
  end
  global _naming_priority_order = np
end

export set_naming_exceptions!

function set_naming_exceptions!(b::Bool)
  _naming_allow_exceptions = b
end

function set_naming_exceptions(d::Dict{Vector{Int64},String})
  !all( c -> length(c) == 4, collect( keys(d)) ) && error("The keys of the dictionary of naming exceptions should be anyonwiki_codes, i.e. integer vectors of lenght 4.")

  global _naming_exceptions = d
end

function name(fr::FusionRing)
  c = anyonwiki_code(fr)
  haskey(_naming_exceptions,c) && return _naming_exceptions[c]

  nms = names(fr)

  for k in _naming_priority_order
    l = nms[k]
    !ismissing(l) && !isempty(l) && return first(l)
  end

  return missing
end

function Base.show(io::IO, ring::FusionRing)
  p(str) = print(io, str);
  n = name(ring)

  !ismissing(n) && return p("FR(" * n * ")")

  if !ismissing(ring.anyonwiki_code)
    return p("FR(" * join( string.(ring.anyonwiki_code) , ", " ) * ")")
  end

  props = string.(comap([rank, multiplicity, nnsd], ring))

  return p("FR(" * join(props, ", ") * ", ? )")
end

export print_multiplication_table

# TO DO: 
# * no nlonger use names - should be labels
# * labels should be printed as bold integers
# """
# Pretty prints the multiplication table as strings (no mutation).
# """
# function print_multiplication_table(fr::FusionRing; include_zeros::Bool=false)
#     N = multiplication_table(fr)
#     names = labels(fr)
#     r = length(names)
#     head = "× │ " * join(names, " │ ")
#     sep  = "──┼" * "───┼"^(r-1) * "──"
#     println(head); println(sep)
#     for i in 1:r
#         rowcells = String[]
#         for j in 1:r
#             d = fusion_product(fr, i, j)
#             if include_zeros
#                 parts = String[]
#                 for c in 1:r
#                     m = get(d,c,0)
#                     if m==0; push!(parts, "0 "*names[c])
#                     elseif m==1; push!(parts, names[c])
#                     else; push!(parts, string(m," ",names[c]))
#                     end
#                 end
#                 push!(rowcells, join(parts, " + "))
#             else
#                 isempty(d) && push!(rowcells, "0") && continue
#                 push!(rowcells,
#                     join([ m==1 ? names[c] : string(m," ",names[c]) for (c,m) in d ], " + "))
#             end
#         end
#         println(names[i], " │ ", join(rowcells, " │ "))
#     end
#     nothing
# end

function print_multiplication_table(r::FusionRing)
  rk = rank(r)
  mt = multiplication_table(r)

  tab = fill("", rk, rk)
  for i in 1:rk, j in 1:rk
    tab[i, j] = row_to_string(r, mt[i, j, :])
  end
  return tab
end

export pmt

pmt = print_multiplication_table

export row_to_string

function row_to_string(r::FusionRing, row)::String
  n            = length(row)
  el_names     = labels(r)
  non_zero_ind = findall(i -> row[i] > 0, 1:n)
  to_string(i) = element_to_string(row[i], el_names[i])

  return join(map(to_string, non_zero_ind), " ⊕ ")
end

function element_to_string(mult, elem)::String
  if mult == 0
    return ""
  elseif mult == 1
    return elem
  else
    return string(mult) * " " * elem
  end
end

"Pretty one-liner: `a × b = ...` using printed names; `a,b` are indices."
function product_string(fr::FusionRing, a::Int, b::Int)
  names = labels(fr)
  d = fusion_product(fr, a, b)
  rhs =
    if isempty(d)
      "0"
    else
      function tostr(tuple)
        c, m = tuple
        m == 1 ? names[c] : string( m, " ", names[c] )
      end

      outcomes = sort(collect(d),by=(x->x[1]))
      join(tostr.(outcomes), " ⊕ ")
    end
  return string(names[a], " × ", names[b], " = ", rhs)
end

function export_tex_reps(filename::String, v::Vector{QQBarFieldElem}; try_cyclo = false)
  data = Dict(qqb_id(x) => tex_reps(x) for x in v)

  open(filename, "w") do f
    return JSON.json(f, data; pretty = true, inline_limit = 10)
  end
end

function tex_reps(x::QQBarFieldElem; try_cyclo = false)
  rat = rational_tex_rep(x)
  if rat != ""
    return Dict(
      "rational"  => rat,
      "radical"   => rat,
      "power_sum" => rat,
      "cyclo"     => rat,
      "general"   => general_tex_rep(x),
    )
  end

  ps = power_sum_tex_rep(x)
  if ps != ""
    cyc = ps
  elseif (ps == "") && try_cyclo
    cyc = cyclo_tex_rep(x)
  else
    cyc = ""
  end

  return Dict(
    "rational"  => rational_tex_rep(x),
    "radical"   => radicals_tex_rep(x),
    "power_sum" => ps,
    "cyclo"     => cyc,
    "general"   => general_tex_rep(x),
  )
end

function general_tex_rep(x::QQBarFieldElem)
  n = string(rootnum(x))

  return fix_poly_string("[" * string(minpoly(x)) * "]_{" * n * "}")
end

function fix_fractions(str::String)::String
  return replace(str, r"(?<p>\d+)//(?<q>\d+)" => s"\\frac{\g<p>}{\g<q>}")
end

function fix_mult(str::String)::String
  return replace(str, "*" => "")
end

function fix_spaces(str::String)::String
  return replace(str, " " => "")
end

function fix_powers(str::String)::String
  return replace(str, r"\^(?<p>\d{2,})" => s"^{\g<p>}")
end

function fix_cyclo(str::String)::String
  return replace(str, r"zeta\((?<deg>\d+)\)" => s"\\zeta_{\g<deg>}")
end

function fix_poly_string(str::String)::String
  return (fix_fractions ∘ fix_mult ∘ fix_spaces ∘ fix_powers)(str)
end

function cyclo_tex_rep(x::QQBarFieldElem)
  cx, emb = to_composite_field(x)
  if !is_abelian(parent(cx))
    return ""
  else
    el = to_cyclotomic_field(cx, emb)[1]

    (fix_cyclo ∘ fix_poly_string ∘ string ∘ QQab)(el)
  end
end

function radicals_tex_rep(x::QQBarFieldElem)
  mp = minpoly(x)
  d  = degree(mp)
  if d == 1
    mp, q = collect(coefficients(mp))
    fix_poly_string(string(-mp//q))
  elseif d == 2
    # output string consists of 4 parts:
    # <-b/2a> * <sign> * <factor/a2> * <sqrt{Δ/factor^2}>
    #   s1    *   s2   *     s3      *        s4
    QQb     = algebraic_closure(QQ)
    c, b, a = collect(coefficients(mp))
    s1      = b == 0 ? "" : string(-b//(2*a))

    Δ    = b^2 - 4*a*c
    pval = (-QQb(b) + sqrt(QQb(Δ)))//(2*a)
    s2   = pval == x ? "+" : "-"

    (d, e) = factor_squares(Δ)
    s3     = d//(2*a) == 1 ? "" : string(d//(2*a))

    s4 = string("\\sqrt{", (Δ > 0 ? e : -e), "}")

    return fix_poly_string(string(s1, s2, s3, s4))
  else
    ""
  end
end

# writes an integer as a multiplication of a squared number
# and a non-square
function factor_squares(x)
  zzx = ZZ(x)
  sq_factors = filter(λ -> λ[2] > 1, (collect ∘ factor)(zzx))
  length(sq_factors) == 0 && return (1, x)

  sq  = prod(a^ZZ(floor(ZZ(b)//2)) for (a, b) in sq_factors)
  nsq = zzx/sq^2

  return (sq, nsq)
end

#TODO: surely there must be a better way...
function rational_tex_rep(x::QQBarFieldElem)
  if isinteger(x)
    m = - first(collect(coefficients(minpoly(x))))
    return string(m)
  elseif is_rational(x)
    cf = collect(coefficients(minpoly(x)))
    return fix_poly_string(string(-cf[1]//cf[2]))
  else
    return ""
  end
end

# represent root of pol of form 1 + ax + (ax)^2 +... +(ax)^n
function power_sum_tex_rep(x::QQBarFieldElem)
  mp = minimal_polynomial(x)
  if is_power_sum(mp)
    gen = (collect(coefficients(mp)))[2]
    n = degree(mp)
    a = x*gen
    i = findfirst(j -> QQb(ζ(n+1)^j) == a, 1:n)
    factorstring = gen == 1 ? "" : string(1//gen)
    return string(factorstring, "\\zeta_{", n+1, "}^{", i, "}")
  else
    return ""
  end
end

function is_power_sum(poly::QQPolyRingElem)
  return (is_geometric_array ∘ collect ∘ coefficients)(poly)
end

function is_geometric_array(a::Vector{T}) where {T}
  if length(a) == 0
    return true
  elseif length(a) == 1
    return a[1] == 1
  else
    gen = a[2]
    for i in 2:(length(a) - 1)
      if gen^i != a[i + 1]
        return false
      end
      continue
    end
    return true
  end
end
