## Needs Hecke 0.39.13

###############################################################################
#
#  Create labels
#
###############################################################################

@doc raw"""
    get_label_and_scaling_factor(G::ZZGenus)

Return the label of the reduced genus of ``G`` and the associated scaling
factor.
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

### Types

@doc raw"""
    DBKeys{T} where T

Object to store an index for the keys of each entry in a database.
"""
mutable struct DBKeys{T}
  _index::Set{T}

  function DBKeys{T}(_index::Set{T}) where T
    return new{T}(_index)
  end
end

@doc raw"""
    ZZLatDefDB

Object to work with the database of definite genera. It stores an absolute
path where to reach the database and an index of keys for the corresponding
dataset.
"""
mutable struct ZZLatDefDB
  path::String
  keys::DBKeys{String}

  function ZZLatDefDB(path)
    z = new(path)
    generate_keys!(z)
    return z
  end
end

### Constructors

@doc raw"""
    definite_lattice_database(path::String) -> ZZLatDefDB

Initiate the object to reach the database of definite genera stored at
the absolute path `path`. If `path` is not set in output, then it chooses by
default the absolute path of the directory in which the file containing this
function is stored.
"""
function definite_lattice_database(path::String = joinpath(@__DIR__, "definite_genera"))
  # Initiate all the necessary folders if not yet existing
  if !isdir(joinpath(path, "dataset"))
    mkdir(joinpath(path, "dataset"))
  end
  if !isdir(joinpath(path, "buffer"))
    mkdir(joinpath(path, "buffer"))
  end
  if !isdir(joinpath(path, "quarantine"))
    mkdir(joinpath(path, "quarantine"))
  end
  if !isdir(joinpath(path, "duplicate"))
    mkdir(joinpath(path, "duplicate"))
  end
  if !isdir(joinpath(path, "temporary"))
    mkdir(joinpath(path, "temporary"))
  end
  return ZZLatDefDB(path)
end

@doc raw"""
    generate_keys!(db::ZZLatDefDB) -> Nothing

Generate the index of keys for the main `dataset` folder of the database of
definite lattices reached by `db` and stores it.
"""
function generate_keys!(db::ZZLatDefDB)
  _path = joinpath(path(db), "dataset")
  _index = Set{String}()
  for label in readdir(_path)
    push!(_index, label)
  end
  ind = DBKeys{String}(_index)
  db.keys = ind
  return nothing
end

### Accessors

@doc raw"""
    path(db::ZZLatDefDB) -> String

Return the absolute path at which the database is stored.
"""
path(db::ZZLatDefDB) = db.path

@doc raw"""
    keys(db::ZZLatDefDB) -> DBKeys

Return the object containing the keys of the entries in the `dataset` folder
in the database reached by `db`.
"""
keys(db::ZZLatDefDB) = db.keys

@doc raw"""
    index(ind::DBKeys{T}) -> Set{T}

Return the index of the set of keys `ind`.
"""
index(ind::DBKeys) = ind._index

@doc raw"""
    index(db::ZZLatDefDB) -> Set{String}

Return the index for the entries in the `dataset` folder in the database
of definite lattices reached by `db`.
"""
index(db::ZZLatDefDB) = index(keys(db))

### Get and set index

# Load the representatives for the genus with symbol G
# Return an error if the genus is not available in the main dataset,
# or if the corresponding entry is corrupted
Base.getindex(db::ZZLatDefDB, G::ZZGenus) = load_genus(db, G)

# Load the entry with key label
# Return an error if the main dataset does not contain such entry,
# or if the entry is corrupted
Base.getindex(db::ZZLatDefDB, label::String) = load_genus(db, label)

# Add a new entry to the main dataset, whose key is the label of the
# reduced genus associated to G and the values are the corresponding
# rescaling of the lattices in lats.
function Base.setindex!(db::ZZLatDefDB, lats::Vector{ZZLat}, G::ZZGenus)
  save_genus!(db, lats, G)
  return nothing
end

### Key handling

@doc raw"""
    add_key!(ind::DBKeys{T}, key::T) -> Nothing

Add `key` to `ind`.
"""
function add_key!(ind::DBKeys{T}, key::T) where T
  push!(index(ind), key)
end

@doc raw"""
    add_key!(db::ZZLatDefDB, label::String) -> Nothing

Add `label` to the keys of the database `db`.
"""
function add_key!(db::ZZLatDefDB, label::String)
  add_key!(keys(db), label)
end

### Containement

# Return whether the label of the reduced genus of `G` is in the
# set of keys encoded in `ind`
function Base.in(G::ZZGenus, ind::DBKeys{String})
  label, _ = get_label_and_scaling_factor(G)
  return label in index(ind)
end

Base.in(G::ZZGenus, db::ZZLatDefDB) = G in keys(db)

# Return whether the reduced genus of `G` is available in the
# database `db`
function Base.haskey(db::ZZLatDefDB, G::ZZGenus)
  return G in keys(db)
end

# Return whether the index of `db` references `label`
function Base.haskey(db::ZZLatDefDB, label::String)
  return label in index(db)
end

### Update function

# Alias for `generate_keys!`
update_keys!(db::ZZLatDefDB) = generate_keys!(db)

###############################################################################
#
#  Safeguards
#
###############################################################################

### Mass exception
# Error if mass formula is not correct

struct MassFormulaError <: Exception
  path_to_genus::String

  function MassFormulaError(path_to_genus::String)
    return new(path_to_genus)
  end
end

function Base.showerror(io::IO, e::MassFormulaError)
  println(io, "!! Database corrupted !! Mass formula is wrong: $(e.path_to_genus) !! Run a cleanup of the main dataset")
end

### Lattice exception
# Error if a lattice file is corrupted

struct LatticeFileError <: Exception
  path_to_file::String
  msg::String

  function LatticeFileError(path_to_file::String, msg::String)
    return new(path_to_file, msg)
  end
end

function Base.showerror(io::IO, e::LatticeFileError)
  println(io, "!! Database corrupted !! $(e.msg): $(e.path_to_file) !! Run a cleanup of the main dataset")
end

### Genus exception
# Error if a genus entry is corrupted

struct GenusDirError <: Exception
  path_to_genus::String
  msg::String

  function GenusDirError(path_to_genus::String, msg::String)
    return new(path_to_genus, msg)
  end
end

function Base.showerror(io::IO, e::GenusDirError)
  println(io, "!! Database corrupted !! $(e.msg): $(e.path_to_genus) !! Run a cleanup of the main dataset")
end

### Checks for corrupted data
# Checks to decide if a certain entry of the database is corrupted
# This should be fast so we check only trivial things

function check_mass_formula(lats::Vector{ZZLat})
  return mass(first(lats)) == sum(Base.Fix1(QQ, 1)∘automorphism_group_order, lats)
end

function is_corrupted_entry(
  db::ZZLatDefDB,
  entry_label::String;
  folder_name::String="dataset",
)
  lats = try
           load_genus(db, entry_label; folder_name)
         catch
           return true
         end

  G = genus(first(lats))
  !all(isequal(G)∘genus, lats) && return true

  _label, _ = get_label_and_scaling_factor(G)
  !contains(entry_label, _label) && return true

  return false
end

### Data validity
# Tests for checking validity of data
# This can be slow, so we add a bit more to the tests

@doc raw"""
    is_valid_entry(
      db::ZZLatDefDB,
      entry_label::String;
      folder_name::String="dataset",
      quarantine::Bool=false,
    )

Decide whether the entry `entry_label` from the folder `folder_name` in the
database of definite lattices reached by `db` is valid.

Here, by valid we mean that the data stored in the corresponding entry defines
a complete set of pairwise non-isometric lattices of the genus with label
`entry_label` which satisfy the mass formula.

If the entry is not valid and `quarantine` is set to `true`, then the entry is
moved to the `quarantine` folder.
"""
function is_valid_entry(
  db::ZZLatDefDB,
  entry_label::String;
  folder_name::String="dataset",
  quarantine::Bool=false,
)
  lats = try
           load_genus(db, entry_label; folder_name)
         catch
           quarantine && move_to_quarantine!(db, folder_name, entry_label)
           return false
         end
  # At that point, all lattices are in the same genus, they know the
  # order of their automorphism group and the mass formula is satisfied
  # It just remains to check that they are really pairwise
  # non-isometric
  for i in 1:length(lats)
    if any(j -> is_isometric(lats[i], lats[j]), i+1:length(lats))
      quarantine && move_to_quarantine!(db, folder_name, entry_label)
      return false
    end
  end
  return true
end

@doc raw"""
    test_validity_database(
      db::ZZLatDefDB,
      test_all::Bool=true,
      quarantine::Bool=true,
    ) -> Int

Test whether all the entries from the main `dataset` folder of the database
of definite lattices reached by `db` are valid. See also `is_valid_entry`.

If an entry is not valid and `quarantine` is set to `true`, then the entry is
moved to the `quarantine` folder.

If `test_all` is set to false, then the function throws an error as soon as
it finds a non-valid entry. Otherwise the function returns the number of
non-valid entries.
"""
function test_validity_database(
  db::ZZLatDefDB;
  test_all::Bool=true,
  quarantine::Bool=false,
)
  update_keys!(db)
  count_nonvalid_entries = Int(0)
  for entry_label in index(db)
    flag = is_valid_entry(db, entry_label; quarantine)
    if !flag
      !test_all && error("Found a corrupted entry")
      count_nonvalid_entries += 1
    end
  end
  return count_nonvalid_entries
end

###############################################################################
#
#  Database management
#
###############################################################################

@doc raw"""
  move_to_quarantine!(
    db::ZZLatDefDB,
    folder_name::String,
    entry_label::String,
  ) -> String

Move the entry `entry_label` from the folder `folder_name` of the database
of definite lattices reached by `db` to the `quarantine` folder.
"""
function move_to_quarantine!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String,
)
  if contains(entry_label, "__")
    k =  findfirst(k -> entry_label[k:k+1] == "__", 1:length(entry_label))
  else
    k = 0
  end
  if !iszero(k)
    trimmed_label = entry_label[1:k-1]
  else
    trimmed_label = entry_label
  end
  _source = joinpath(path(db), folder_name, entry_label)
  @assert isdir(_source)
  j = count(contains(trimmed_label), readdir(joinpath(path(db), "quarantine")))
  quarantine_label = trimmed_label*"__qua.$(j)"
  _dest = joinpath(path(db), "quarantine", quarantine_label)
  @assert !isdir(_dest)
  mv(_source, _dest)
  return quarantine_label
end

@doc raw"""
    cleanup_database!(db::ZZLatDefDB) -> Nothing

Go through the main folder `dataset` of the database of definite lattices
reached by `db` and put in quarantine any corrupted entry.
"""
function cleanup_database!(db::ZZLatDefDB)
  update_keys!(db)
  for label in index(db)
    if is_corrupted_entry(db, label)
      move_to_quarantine!(db, "dataset", label)
    end
  end
  return nothing
end

@doc raw"""
    empty_folder!(db::ZZLatDefDB, folder_name::String) -> Nothing

Remove all entries the folder `folder_name` in the database of definite
lattices reached by `db`.

For data safety, this is only available for all folders except the main folder
`dataset`.
"""
function empty_folder!(db::ZZLatDefDB, folder_name::String)
  @assert folder_name != "dataset"
  _path = joinpath(path(db), folder_name)
  !isdir(_path) && return nothing
  for entry_label in readdir(_path)
    delete_from_folder!(db, folder_name, entry_label)
  end
  return nothing
end

@doc raw"""
    delete_from_folder!(
      db::ZZLatDefDB,
      folder_name::String,
      entry_label::String
    ) -> Nothing

Delete the entry `entry_label` in the folder `folder_name` in the database
of definite lattices reached by `db`.

For data safety, this is only available for all folders except the main folder
`dataset`.
"""
function delete_from_folder!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String,
)
  @assert folder_name != "dataset"
  _path = joinpath(path(db), folder_name)
  !isdir(_path) && return nothing
  _genus = joinpath(_path, entry_label)
  if isdir(_genus)
    rm(_genus; recursive=true)
  end
  return nothing
end

@doc raw"""
    promote_entry_to_main_dataset!(
      db::ZZLatDefDB,
      folder_name::String,
      entry_label::String;
      replace_old::Bool=false,
      ignore_conflict::Bool=true,
      quarantine::Bool=true,
    ) -> Nothing

Move the entry `entry_label` from the folder `folder_name` in the database of
definite lattices reached by `db` to the main folder `dataset`.

If the entry is already available in `dataset`, then
  - either `replace_old` is set to true and the content is replaced by the one
    of `entry_label` from `folder_name`;
  - or `ignore_conflict` is set to `true` and the content of `entry_label` from
    `folder_name` is deleted;
  - or an error is thrown.

This function cannot be called for context in quarantine. If the content of
`entry_label` from `folder_name` is corrupted, then:
  - either `quarantine` is set to true and the data is moved in the
    `quarantine` folder,
  - or an error is thrown.
"""
function promote_entry_to_main_dataset!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String;
  replace_old::Bool=false,
  ignore_conflict::Bool=true,
  quarantine::Bool=true,
)
  @assert folder_name in ["buffer", "duplicate"]
  if is_corrupted_entry(db, entry_label; folder_name)
    if quarantine
      move_to_quarantine!(db, folder_name, entry_label)
      return nothing
    else
      error("Corrupted entry: "*entry_label)
    end
  end

  _source = joinpath(path(db), folder_name, entry_label)
  @assert isdir(_source)
  k = findfirst(k -> entry_label[k:k+1] == "__", 1:length(entry_label))
  @assert !isnothing(k)
  trimmed_label = entry_label[1:k-1]
  _dest = joinpath(path(db), "dataset", trimmed_label)

  if isdir(_dest)
    if replace_old
      mv(_source, _dest; force=true)
      return nothing
    elseif ignore_conflict
      rm(_source; recursive=true)
      return nothing
    else
      error("Data already available in the main dataset")
    end
  end

  mv(_source, _dest)
  add_key!(db, entry_label)
  return nothing
end

@doc raw"""
    promote_folder_to_main_dataset!(
      db::ZZLatDefDB,
      folder_name::String,
      replace_old::Bool=false,
      ignore_conflict::Bool=true,
      quarantine::Bool=true,
    ) -> Nothing

Move all the entries from the folder `folder_name` in the database of
definite lattices reached by `db` to the main folder `dataset`.

See `promote_entry_to_main_dataset!` for further documentation.
"""
function promote_folder_to_main_dataset!(
  db::ZZLatDefDB,
  folder_name::String;
  replace_old::Bool=false,
  ignore_conflict::Bool=true,
  quarantine::Bool=true,
)
  _path = joinpath(path(db), folder_name)
  !isdir(_path) && return nothing
  for entry_label in readdir(_path)
    promote_entry_to_main_dataset!(db, folder_name, entry_label; replace_old, ignore_conflict, quarantine)
  end
  return nothing
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
  length(V) == binomial(n+1, 2) || throw(LatticeFileError(lat, "Incompatible rank and half gram"))
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

If ``G`` can be loaded from the main `dataset` folder of the database of
definite lattices reached by `db`, return the corresponding representatives for
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
      entry_label::String,
      s::QQFieldElem = QQ(1);
      folder_name::String="dataset",
    ) -> Vector{ZZLat}

Return the list of lattices ``L(s)`` where ``L`` ranges over all the lattices
stored at the entry `entry_name` from the folder `folder_name` in the database
of definite lattices reached by `db`.

The keyword argument `folder_name` refers to the main `dataset` folder by
default, but one can also input any other folder which is not `quarantine`.

If the entry `entry_label` from the folder `folder_name` does not exist, is
empty or is corrupted, an error is thrown.
"""
function load_genus(
  db::ZZLatDefDB,
  entry_label::String,
  s::QQFieldElem = QQ(1);
  folder_name::String="dataset",
)
  lats = ZZLat[]
  _path = joinpath(path(db), folder_name, entry_label)
  @req isdir(_path) "Entry does not exist in the given folder"
  rd = readdir(_path; join=true)
  isempty(rd) && throw(GenusDirError(_path, "Empty genus entry"))
  for lat in rd
    contains(lat, "lat_") && endswith(lat, ".txt") || throw(GenusDirError(_path, "Wrong lattice file format"))
    push!(lats, load_lattice(lat, s))
  end

  !check_mass_formula(lats) && throw(MassFormulaError(_path))
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
  length(data) < 3 && throw(LatticeFileError(lat, "Missing entries in lattice file"))
  rk = try
         Base.parse(Int, first(data))
       catch
         throw(LatticeFileError(lat, "Cannot parse first line"))
       end

  _, V = try
           Hecke._parse(Vector{QQFieldElem}, IOBuffer(data[2]))
         catch
           throw(LatticeFileError(lat, "Cannot parse second line"))
         end

  L = lattice_from_data(V, rk, s)

  L.automorphism_group_order = try
                                 last(Hecke._parse(ZZRingElem, IOBuffer(data[3])))
                               catch
                                 throw(LatticeFileError(lat, "Cannot parse third line"))
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
      s::QQFieldElem = QQ(1);
    ) -> String

Given an integer lattice ``L = M(s)`` where ``M`` is indivisible integral of
nonnegative signature, return a string to be written in a ".txt" file to
store the lattice ``M``. The string describes 3 lines:
- the first line encodes the rank of ``M``,
- the second line encodes the entries of the upper triangular part of the
  Gram matrix of ``M``,
- the third line encodes the order of the automorphism group of ``L``
  (which is the same as for ``M``).
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
  s = automorphism_group_order(L)
  str *= "\n$(s)"
  return str
end

@doc raw"""
    save_genus!(
      db::ZZLatDefDB,
      lats::Vector{ZZLat},
      G::ZZGenus = genus(first(lats));
      safe::Bool=true,
      replace_corrupted_data::Bool=true,
      save_duplicate::Bool=false,
      verbose::Bool=false,
    ) -> Nothing

Update the main `dataset` folder of the database of definite lattice reached by
`db` by storing the (rescaled) content of `lats` in a new entry labeled by
the reduced of ``G``. Other ways to call this function:
- `setindex!(db, lats, G)`
- `db[G] = lats`

If `safe` is to `true`, the function checks whether the data is safe to be
saved (highly recommended). It essentially checks:
  - whether all the lattices in `lats` are of the genus `G`,
  - whether the mass formula for `G` is satisfied by the lattices in `lats`.
If the data is not safe, it is stored in the `quarantine` folder.

If the genus ``G`` is already available in `dataset`, then:
  - either the corresponding entry in `dataset` is corrupted. In this case, it
    is moved to `quarantine` and
      * either `replace_corrupted_data` is `true` and the new entry is stored
        in the main `dataset` folder.
      * or the new entry is stored in the folder `buffer`.
  - or `save_duplicate` is set to `true` and the new entry is saved in the
    folder `duplicate`.
  - or the new entry is ignored.

If `verbose` is set to `true`, the user allows the function to tell which
decision is taken for storing the new entry.
"""
function save_genus!(
  db::ZZLatDefDB,
  lats::Vector{ZZLat},
  G::ZZGenus = genus(first(lats));
  safe::Bool=true,
  replace_corrupted_data::Bool=true,
  save_duplicates::Bool=false,
  verbose::Bool=false,
)
  label, s = get_label_and_scaling_factor(G)
  verbose && println("Store new genus: ", label)
  add_new_key = true

  if safe
    if !all(isequal(G)∘genus, lats) || !check_mass_formula(lats)
      add_new_key = false
      verbose && println("!! Data corrupted !! The lattices are not all in the same genus or do not satisfy the mass formula...")
      j = count(contains(label), readdir(joinpath(path(db), "quarantine")))
      quarantine_label = label*"__qua.$(j)"
      _path = joinpath(path(db), "quarantine", quarantine_label)
      verbose && println(">>>> New data stored in quarantine folder: ", quarantine_label)
    end
  end

  if add_new_key
    _path = joinpath(path(db), "dataset", label)
    if isdir(_path) && !save_duplicates && !is_corrupted_entry(db, label)
      verbose && println(">>>> New data is ignored")
      return true
    end
  end

  tmp_path = mktempdir(joinpath(path(db), "temporary"); prefix=label)
  for i in 1:length(lats)
    lat = joinpath(tmp_path, "lat_$(i).txt")
    touch(lat)
    _f = open(lat, "w")
    Base.write(_f, data_from_lattice(lats[i], s))
    close(_f)
  end

  if is_corrupted_entry(db, tmp_path; folder_name="temporary")
    add_new_key = false
    verbose && println("!! Data corrupted !! Something went wrong during saving...")
    j = count(contains(label), readdir(joinpath(path(db), "quarantine")))
    quarantine_label = label*"__qua.$(j)"
    _path = joinpath(path(db), "quarantine", quarantine_label)
    verbose && println(">>>> New data stored in quarantine folder: ", quarantine_label)
  end

  if isdir(_path)
    verbose && println("!! Duplicate warning !! The genus seems to be already in the database...")
    if is_corrupted_entry(db, label)
      verbose && println(">> Current version of genus is corrupted...")
      quarantine_label = move_to_quarantine!(db, "dataset", label)
      verbose && println(">>>> Old data moved to quarantine folder: ", quarantine_label)
      if !replace_corrupted_data
        add_new_key = false
        j = count(contains(label), readdir(joinpath(path(db), "buffer")))
        buffer_label = label*"__buf.$(j)"
        _path = joinpath(path(db), "buffer", buffer_label)
        verbose && println(">>>> New data stored in buffer folder: ", buffer_label)
      else
        verbose && println(">>>> New data stored in dataset folder: ", label)
      end
    else
      verbose && println(">> Current version of genus is not corrupted...")
      if save_duplicates
        j = count(contains(label), readdir(joinpath(path(db), "duplicate")))
        duplicate_label = label*"__dup.$(j)"
        _path = joinpath(path(db), "duplicate", duplicate_label)
        verbose && println(">>>> New data stored in duplicate folder: ", duplicate_label)
      else
         verbose && println(">>>> New data is ignored")
         return true
      end
    end
  elseif add_new_key
    verbose && println(">>>> New data stored in dataset folder: ", label)
  end

  mv(tmp_path, _path)
  add_new_key && add_key!(db, label)
  return add_new_key
end

###############################################################################
#
#  Temporary script: move old database to new infrastructure
#
#  To be removed once the old data has been transfered in
#  `definite_genera/dataset`
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

function move_to_new_database(
  db_path::String = @__DIR__;
  verbose::Bool=false,
)
  db = definite_lattice_database(joinpath(db_path, "definite_genera"))
  rd = readdir(joinpath(db_path, "defgen_db"); join=true)
  for rank_folder in rd
    for genus_entry in readdir(rank_folder; join=true)
      lats = load_genus(genus_entry)
      G = genus(first(lats))
      if haskey(db, G)
        continue
      end
      flag = save_genus!(db, lats, G; verbose)
      @assert !flag || Set(gram_matrix.(lats)) == Set(gram_matrix.(db[G]))
    end
  end
  return nothing
end
