###############################################################################
#
#  Compact storage of a genus
#
###############################################################################

# M must be symmetric
function _get_half_gram(L::ZZLat)
  M = gram_matrix(L)
  str = "$(nrows(M))\n["
  for i in 1:nrows(M), j in i:ncols(M)
    str *= "$(M[i,j]),"
  end
  str = str[1:end-1]*"]"
  if isdefined(L, :automorphism_group_order)
    str *= "\n$(L.automorphism_group_order))"
  end
  return str
end

# V is a list of gram matrices of the lattices in the genus
# f is a path to a numbered dir in the correct rank folder
function save_genus(V::Vector{ZZLat}, f::String)
  for i in 1:length(V)
    p = f*"/lat_$(i).txt"
    touch(p)
    _f = open(p, "w")
    Base.write(_f, _get_half_gram(V[i]))
    close(_f)
  end
  return nothing
end

###############################################################################
#
#  Reading a stored genus
#
###############################################################################

# V is a list of integers, half a gram matrix
# n is the rank of the matrix to be constructed
function _gram_from_list(V::Vector{QQFieldElem}, n::Int)
  M = zero_matrix(QQ, n, n)
  k = 0
  for i in 1:n
    for j in i:n
      k += 1
      M[i, j] = M[j, i] = V[k]
    end
  end
  return integer_lattice(; gram=M)
end

# Load the genus stored in the numbered dir f
function load_genus(f::String)
  gg = ZZLat[]
  files = readdir(f; join=true)
  for file in files
    rl = readlines(file)
    n = Base.parse(Int, rl[1])
    _, V = Hecke._parse(Vector{QQFieldElem}, IOBuffer(rl[2]))
    L = _gram_from_list(V, n)
    if length(rl) == 3
      L.automorphism_group_order = Hecke._parse(ZZRingElem, IOBuffer(rl[3]))[2]
    end
    push!(gg, L)
  end
  return gg
end

###############################################################################
#
#  Load the database of definite genera
#
###############################################################################

# The output is a dictionary ZZGenus -> list of gram matrices
function load_db(f::String)
  D = Dict{ZZGenus, Vector{ZZLat}}()
  df = readdir(f; join=true)
  for _p in df
    rf = readdir(_p; join=true)
    for __p in rf
      gg = load_genus(__p)
      is_empty(gg) && continue
      G = genus(gg[1])
      @assert !haskey(D, G)
      D[G] = gg
    end
  end
  return D
end

###############################################################################
#
#  Update database
#
##############################################################################

function update_db!(f::String, D::Dict{ZZGenus, Vector{ZZLat}})
  D2 = load_db(f)
  c = collect(keys(D))
  c2 = collect(keys(D2))
  setdiff!(c, c2)
  
  for G in c
    add_to_db!(f, D[G])
  end
  return nothing
end

function add_to_db!(f::String, V::Vector{ZZLat})
  @assert !isempty(V)
  G = genus(V[1])
  r = rank(G)
  fr = f*"rank$r/"
  if !isdir(fr)
    mkdir(fr)
  end
  l = length(readdir(fr))
  n = 5-ndigits(l+1)
  frl = fr*"0"^n*"$(l+1)"
  mkdir(frl)
  save_genus(V, frl)
  return nothing
end


