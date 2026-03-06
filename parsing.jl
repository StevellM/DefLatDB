## Needs Hecke 0.39.13

###############################################################################
#
#  Create labels
#
###############################################################################

@doc raw"""
    get_label_and_scaling_factor(G::ZZGenus)

Return the label of the reduced of ``G`` and the associated scaling factor.
"""
function get_label_and_scaling_factor(G::ZZGenus)
  Gred, scaling_factor = Hecke.reduced_genus_with_data(G)
  t, r, s, d, l, symbolGred = Hecke.genus_info(Gred)
  label = "$r.$s.$d.$l.$t"*symbolGred
  return label, scaling_factor
end

###############################################################################
#
#  Database object
#
###############################################################################

@doc raw"""
    ZZLatDefDBKeys

Type to store unique keys for reduced genera in the database of definite
genera.
"""
struct ZZLatDefDBKeys
  keys::Set{String}

  function ZZLatDefDBKeys(keys::Set{String})
    return new(keys)
  end
end

@doc raw"""
    ZZLatDefDB

Type to work with database of definite genera.
"""
struct ZZLatDefDB
  path::String
  keys::ZZLatDefDBKeys

  function ZZLatDefDB(path)
    keys = generate_keys(path)
    return new(path, keys)
  end
end

@doc raw"""
    definite_lattice_database(path::String) -> ZZLatDefDB

Initiate the database object for the database of definite genera stored at
the absolute path `path`. If `path` is not set in output, then it chooses by
default the absolute path of the directory in which the file containing this
function is stored.
"""
function definite_lattice_database(path::String = joinpath(@__DIR__, "definite_genera"))
  return ZZLatDefDB(path)
end

@doc raw"""
    path(db::ZZLatDefDB) -> String

Return the absolute path at which the database is stored.
"""
path(db::ZZLatDefDB) = db.path

@doc raw"""
    keys(db::ZZLatDefDB) -> ZZLatDefDBKeys

Return the object containing the keys of unique elements in `db`.
"""
keys(db::ZZLatDefDB) = db.keys

Base.getindex(db::ZZLatDefDB, G::ZZGenus) = load_genus(db, G)

Base.setindex!(db::ZZLatDefDB, lats::Vector{ZZLat}, G::ZZGenus) = save_genus!(db, lats, G)

function generate_keys(path::String)
  k = Set{String}()
  for label in readdir(path)
    push!(k, label)
  end
  return ZZLatDefDBKeys(k)
end

function add_key!(k::ZZLatDefDBKeys, label::String)
  push!(k.keys, label)
end

function Base.in(G::ZZGenus, k::ZZLatDefDBKeys)
  label, _ = get_label_and_scaling_factor(G)
  return label in k.keys
end

function Base.haskey(db::ZZLatDefDB, G::ZZGenus)
  return G in keys(db)
end

###############################################################################
#
#  Load a definite genus
#
###############################################################################

@doc raw"""
    lattice_from_data(
      V::Vector{QQFieldElem},
      n::Int,
      s::QQFieldElem,
    ) -> ZZLat

Return the integer lattice ``L(s)`` where ``L`` has rank ``n`` and the upper
triangular part of the Gram matrix of ``L`` is given by the entries in `V`.
"""
function lattice_from_data(
  V::Vector{QQFieldElem},
  n::Int,
  s::QQFieldElem = QQ(1),
)
  if !isone(s)
    map!(Base.Fix2(mul!, s), V)
  end
  M = zero_matrix(QQ, n, n)
  k = 0
  for i in 1:n
    for j in i:n
      k += 1
      M[i, j] = M[j, i] = V[k]
    end
  end
  return integer_lattice(; gram=M, cached=false)
end

@doc raw"""
    load_genus(
      db::ZZLatDefDb,
      G::ZZGenus,
    ) -> Vector{ZZLat}

If ``G`` can be loaded from `db`, return the corresponding representatives for
the isometry classes in ``G``.

Other ways to call this function
- `getindex(db, G)`
- `db[G]`
"""
function load_genus(
  db::ZZLatDefDB,
  G::ZZGenus,
)
  label, s = get_label_and_scaling_factor(G)
  return load_genus(db, label, s)
end

@doc raw"""
    load_genus(
      db::ZZLatDefDB,
      label::String,
      s::QQFieldElem = QQ(1),
    ) -> Vector{ZZLat}

Return the list of lattices ``L(s)`` where ``L`` is one of the lattices
stored in the genus in `db` with label `label`.
"""
function load_genus(
  db::ZZLatDefDB,
  label::String,
  s::QQFieldElem = QQ(1),
)
  lats = ZZLat[]
  _path = path(db)
  _path = joinpath(_path, label)
  @req isdir(_path) "Not in the database"
  rd = readdir(_path; join=true)
  @assert !isempty(rd)
  for lat in rd
    contains(lat, "lat_") && endswith(lat, ".txt") || continue
    push!(lats, load_lattice(lat, s))
  end
  return lats
end

@doc raw"""
    load_lattice(
      lat::String,
      s::QQFieldElem = QQ(1),
    ) -> ZZLat

Return the lattice ``L(s)`` where ``L`` is the lattice stored in the file
with absolute path `lat`.
"""
function load_lattice(
  lat::String,
  s::QQFieldElem = QQ(1),
)
  data = readlines(lat)
  @assert length(data) >= 2
  rk = Base.parse(Int, first(data))
  _, V = Hecke._parse(Vector{QQFieldElem}, IOBuffer(data[2]))
  L = lattice_from_data(V, rk, s)
  if length(data) >= 3
    L.automorphism_group_order = last(Hecke._parse(ZZRingElem, IOBuffer(data[3])))
  end
  return L
end

###############################################################################
#
#  Save a definite genus
#
###############################################################################

@doc raw"""
    data_from_lattice(
      L::ZZLat,
      s::QQFieldElem = QQ(1),
    ) -> String

Given an integer lattice ``L = M(s)`` where ``M`` is indivisible integral of
nonnegative signature, return a string to be written in a ".txt" file to
store the lattice ``M``. The string described 2 or 3 lines:
- the first line encodes the rank of ``M``,
- the second line encodes the entries of the upper triangular part of the
  Gram matrix of ``M``,
- if it has been computed, the third line encodes the order of the
  automorphism group of ``L`` (which is the same as for ``M``).
"""
function data_from_lattice(L::ZZLat, s::QQFieldElem = QQ(1))
  M = gram_matrix(L)
  str = "$(nrows(M))\n["
  z = QQ(0)
  for i in 1:nrows(M), j in i:ncols(M)
    k = M[i,j]
    if !isone(s)
      mul!(z, k, 1//s)
    end
    str *= "$z,"
  end
  str = str[1:end-1]*"]"
  if isdefined(L, :automorphism_group_order)
    str *= "\n$(L.automorphism_group_order)"
  end
  return str
end

@doc raw"""
    save_genus!(
      db::ZZLatDefDB,
      lats::Vector{ZZLat},
      G::ZZGenus = genus(first(lats)),
    ) -> Nothing

Update the database reached by `db` by adding the new genus ``G``, if not
yet available. The associated lattices are given in `lats`.

Other ways to call this function:
- `setindex!(db, lats, G)`
- `db[G] = lats`
"""
function save_genus!(
  db::ZZLatDefDB,
  lats::Vector{ZZLat},
  G::ZZGenus = genus(first(lats)),
)
  @assert all(isequal(G)∘genus, lats)
  add_new_key = true
  label, s = get_label_and_scaling_factor(G)
  _path = joinpath(path(db), label)
  if isdir(_path)
    print("Warning: Definite genus already found in database")
    return nothing
  else
    mkdir(_path)
  end
  for i in 1:length(lats)
    lat = joinpath(_path, "lat_$(i).txt")
    @assert !isfile(lat)
    touch(lat)
    _f = open(lat, "w")
    Base.write(_f, data_from_lattice(lats[i], s))
    close(_f)
  end
  add_new_key && add_key!(keys(db), label)
  return nothing
end

###############################################################################
#
#  Temporary script
#
###############################################################################

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
  return integer_lattice(; gram=M, cached=false)
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

function move_to_new_database(db_path::String = @__DIR__)
  db = definite_lattice_database(joinpath(db_path, "definite_genera"))
  rd = readdir(joinpath(db_path, "defgen_db"); join=true)[end:end]
  for path_to_rank in rd
    for path_to_genus in readdir(path_to_rank; join=true)
      lats = load_genus(path_to_genus)
      G = genus(first(lats))
      if haskey(db, G)
        @assert length(db[G]) == length(lats)
        continue
      end
      db[G] = lats
      @assert Set(gram_matrix.(lats)) == Set(gram_matrix.(db[G]))
    end
  end
  return nothing
end
