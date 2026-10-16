ads###############################################################################
#
#  Create labels
#
###############################################################################

@doc raw"""
    get_label_and_scaling_factor(G::ZZGenus) -> String, QQFieldElem

Return the label of the reduced genus ``G_{red}`` of ``G`` and the scaling
factor ``s`` such that ``G`` is the rescaling ``G_{red}(s)``.

The reduced genus is the rescaling of ``G`` which is integral, primitive (the
scale is ``1``) and of nonnegative signature, so ``s`` is negative for a
negative definite genus.
"""
function get_label_and_scaling_factor(G::ZZGenus)
  Gred, scaling_factor = Hecke.reduced_genus_with_data(G)
  t, r, s, d, l, symbolGred = Hecke.genus_info(Gred)
  label = "$r.$s.$d.$l.$t"*symbolGred
  return label, scaling_factor
end

# An entry outside of the `main` folder is labeled `<label>__<suffix>`
function _trim_label(entry_label::String)
  k = findfirst("__", entry_label)
  return isnothing(k) ? entry_label : entry_label[1:first(k)-1]
end

# Whether the entry `entry_label` stores the genus with label `label`
_is_entry_of(entry_label::String, label::String) = _trim_label(entry_label) == label

# The rank of the genus stored at `entry_label` (the first field of the label)
function _rank_from_label(entry_label::String)
  j = findfirst(isequal('.'), entry_label)
  n = isnothing(j) ? nothing : tryparse(Int, entry_label[1:j-1])
  isnothing(n) && throw(ArgumentError("not a label of the database: $entry_label"))
  return n
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

  function ZZLatDefDB(path::String)
    z = new(path)
    generate_keys!(z)
    return z
  end
end

### Constructors

# the data shipped with this repository
const DEFAULT_DATA_PATH = normpath(joinpath(@__DIR__, "..", "data"))

const FOLDERS = ("main", "quarantine", "duplicate", "temporary")

@doc raw"""
    definite_lattice_database(path::String = DefLatDB.DEFAULT_DATA_PATH) -> ZZLatDefDB

Initiate the object to reach the database of definite genera stored at the
path `path`. By default, this is the folder `data` of this repository.

The folders `main`, `quarantine`, `duplicate` and `temporary` are created in
`path` if they do not exist yet.
"""
function definite_lattice_database(path::String = DEFAULT_DATA_PATH)
  for folder in FOLDERS
    mkpath(joinpath(path, folder))
  end
  return ZZLatDefDB(path)
end

@doc raw"""
    generate_keys!(db::ZZLatDefDB) -> Nothing

Generate the index of keys for the `main` folder of the database of
definite lattices reached by `db` and stores it.
"""
function generate_keys!(db::ZZLatDefDB)
  _path = joinpath(path(db), "main")
  _index = Set{String}(filter(x -> !startswith(x, "."), readdir(_path)))
  db.keys = DBKeys{String}(_index)
  return nothing
end

### Accessors

@doc raw"""
    path(db::ZZLatDefDB) -> String

Return the path at which the database is stored.
"""
path(db::ZZLatDefDB) = db.path

@doc raw"""
    Base.keys(db::ZZLatDefDB) -> DBKeys

Return the object containing the keys of the entries in the `main` folder
in the database reached by `db`.
"""
Base.keys(db::ZZLatDefDB) = db.keys

@doc raw"""
    index(ind::DBKeys{T}) -> Set{T}

Return the index of the set of keys `ind`.
"""
index(ind::DBKeys) = ind._index

@doc raw"""
    index(db::ZZLatDefDB) -> Set{String}

Return the index for the entries in the `main` folder in the database
of definite lattices reached by `db`.
"""
index(db::ZZLatDefDB) = index(keys(db))

Base.length(ind::DBKeys) = length(index(ind))
Base.length(db::ZZLatDefDB) = length(keys(db))

### Get and set index

# Load the representatives for the genus G
# Return an error if the genus is not available in the `main` folder,
# or if the corresponding entry is corrupted
Base.getindex(db::ZZLatDefDB, G::ZZGenus) = load_genus(db, G)

# Load the entry with key label
# Return an error if the `main` folder does not contain such entry,
# or if the entry is corrupted
Base.getindex(db::ZZLatDefDB, label::String) = load_genus(db, label)

# Add a new entry to the `main` folder, whose key is the label of the
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
  return nothing
end

@doc raw"""
    add_key!(db::ZZLatDefDB, label::String) -> Nothing

Add `label` to the keys of the database `db`.
"""
add_key!(db::ZZLatDefDB, label::String) = add_key!(keys(db), label)

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
Base.haskey(db::ZZLatDefDB, G::ZZGenus) = G in keys(db)

# Return whether the index of `db` references `label`
Base.haskey(db::ZZLatDefDB, label::String) = label in index(db)

### Update function

@doc raw"""
    update_keys!(db::ZZLatDefDB) -> Nothing

Alias for `generate_keys!`: read the index of the `main` folder again (e.g.
to see the entries that other processes have stored meanwhile).
"""
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
end

function Base.showerror(io::IO, e::MassFormulaError)
  print(io, "!! Database corrupted !! Mass formula is wrong: $(e.path_to_genus) !! Run a cleanup of the main folder")
end

### Genus exception
# Error if a genus file is corrupted

struct GenusError <: Exception
  path_to_genus::String
  msg::String
end

function Base.showerror(io::IO, e::GenusError)
  print(io, "!! Database corrupted !! $(e.msg): $(e.path_to_genus) !! Run a cleanup of the main folder")
end

### Checks for corrupted data
# Checks to decide if a certain entry of the database is corrupted
# This should be fast so we check only trivial things

function check_mass_formula(lats::Vector{ZZLat})
  return mass(first(lats)) == sum(Base.Fix1(QQ, 1)∘automorphism_group_order, lats; init=QQ(0))
end

@doc raw"""
    is_corrupted_entry(
      db::ZZLatDefDB,
      entry_label::String;
      folder_name::String="main",
    ) -> Bool

Return whether the entry `entry_label` of the folder `folder_name` of the
database `db` cannot be loaded, violates the mass formula, or stores a genus
whose label is not `entry_label` (up to a suffix `__...`). These checks are
cheap; see `is_valid_entry` for a complete check.
"""
function is_corrupted_entry(
  db::ZZLatDefDB,
  entry_label::String;
  folder_name::String="main",
)
  lats = try
           load_genus(db, entry_label; folder_name)
         catch
           return true
         end
  G = genus(first(lats))
  _label, _ = get_label_and_scaling_factor(G)
  return !_is_entry_of(entry_label, _label)
end

### Data validity
# Tests for checking validity of data
# This can be slow, so we add a bit more to the tests

@doc raw"""
    is_valid_entry(
      db::ZZLatDefDB,
      entry_label::String;
      folder_name::String="main",
      quarantine::Bool=false,
    ) -> Bool

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
  folder_name::String="main",
  quarantine::Bool=false,
)
  if is_corrupted_entry(db, entry_label; folder_name)
    quarantine && move_to_quarantine!(db, folder_name, entry_label)
    return false
  end
  # At that point, the entry loads, the lattices know the order of their
  # automorphism group and the mass formula is satisfied. It remains to
  # check that they all lie in the same genus and are pairwise
  # non-isometric.
  lats = load_genus(db, entry_label; folder_name)
  G = genus(first(lats))
  if !all(isequal(G)∘genus, lats)
    quarantine && move_to_quarantine!(db, folder_name, entry_label)
    return false
  end

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
      db::ZZLatDefDB;
      test_all::Bool=true,
      quarantine::Bool=false,
    ) -> Int

Test whether all the entries from the `main` folder of the database
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
  count_nonvalid_entries = 0
  for entry_label in collect(index(db))
    flag = is_valid_entry(db, entry_label; quarantine)
    if !flag
      !test_all && error("Found a corrupted entry: $entry_label")
      count_nonvalid_entries += 1
    end
  end
  quarantine && update_keys!(db)
  return count_nonvalid_entries
end

###############################################################################
#
#  Database management
#
###############################################################################

# A label `<label>__<tag>.<j>` that is not used in `folder_name` yet
function _free_label(db::ZZLatDefDB, folder_name::String, label::String, tag::String)
  j = 0
  while isfile(joinpath(path(db), folder_name, "$(label)__$(tag).$(j)"))
    j += 1
  end
  return "$(label)__$(tag).$(j)"
end

@doc raw"""
    move_to_quarantine!(
      db::ZZLatDefDB,
      folder_name::String,
      entry_label::String,
    ) -> String

Move the entry `entry_label` from the folder `folder_name` of the database
of definite lattices reached by `db` to the `quarantine` folder, and return
its label there.
"""
function move_to_quarantine!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String,
)
  _source = joinpath(path(db), folder_name, entry_label)
  @req isfile(_source) "Entry $entry_label does not exist in the folder $folder_name"
  quarantine_label = _free_label(db, "quarantine", _trim_label(entry_label), "qua")
  mv(_source, joinpath(path(db), "quarantine", quarantine_label))
  folder_name == "main" && delete!(index(db), entry_label)
  return quarantine_label
end

@doc raw"""
    cleanup_database!(db::ZZLatDefDB) -> Nothing

Go through the `main` folder of the database of definite lattices
reached by `db` and put in quarantine any corrupted entry.
"""
function cleanup_database!(db::ZZLatDefDB)
  update_keys!(db)
  for label in collect(index(db))  # collect because the loop changes the index set we loop over
    if is_corrupted_entry(db, label)
      move_to_quarantine!(db, "main", label)
    end
  end
  return nothing
end

@doc raw"""
    empty_folder!(db::ZZLatDefDB, folder_name::String) -> Nothing

Remove all entries the folder `folder_name` in the database of definite
lattices reached by `db`.

For data safety, this is only available for all folders except the `main` folder.
"""
function empty_folder!(db::ZZLatDefDB, folder_name::String)
  @req folder_name != "main" "The main folder cannot be emptied"
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

For data safety, this is only available for all folders except the `main` folder.
"""
function delete_from_folder!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String,
)
  @req folder_name != "main" "Entries of the main folder cannot be deleted"
  _genus = joinpath(path(db), folder_name, entry_label)
  isfile(_genus) && rm(_genus)
  return nothing
end

@doc raw"""
    promote_entry_to_main_folder!(
      db::ZZLatDefDB,
      folder_name::String,
      entry_label::String;
      replace_old::Bool=false,
      ignore_conflict::Bool=true,
      quarantine::Bool=true,
    ) -> Nothing

Move the entry `entry_label` from the folder `folder_name` (currently only
`duplicate`) in the database of definite lattices reached by `db` to the
`main` folder.

If the entry is already available in `main`, then
  - either `replace_old` is set to true and the content is replaced by the one
    of `entry_label` from `folder_name`;
  - or `ignore_conflict` is set to `true` and the content of `entry_label` from
    `folder_name` is deleted;
  - or an error is thrown.

This function cannot be called for content in quarantine. If the content of
`entry_label` from `folder_name` is corrupted, then:
  - either `quarantine` is set to true and the data is moved in the
    `quarantine` folder,
  - or an error is thrown.
"""
function promote_entry_to_main_folder!(
  db::ZZLatDefDB,
  folder_name::String,
  entry_label::String;
  replace_old::Bool=false,
  ignore_conflict::Bool=true,
  quarantine::Bool=true,
)
  @req folder_name == "duplicate" "Only entries of the duplicate folder can be promoted"
  if is_corrupted_entry(db, entry_label; folder_name)
    if quarantine
      move_to_quarantine!(db, folder_name, entry_label)
      return nothing
    else
      error("Corrupted entry: "*entry_label)
    end
  end

  _source = joinpath(path(db), folder_name, entry_label)
  trimmed_label = _trim_label(entry_label)
  _dest = joinpath(path(db), "main", trimmed_label)

  if isfile(_dest)
    if replace_old
      mv(_source, _dest; force=true)
      add_key!(db, trimmed_label)
      return nothing
    elseif ignore_conflict
      rm(_source)
      return nothing
    else
      error("Data already available in the main folder")
    end
  end

  mv(_source, _dest)
  add_key!(db, trimmed_label)
  return nothing
end

@doc raw"""
    promote_folder_to_main_folder!(
      db::ZZLatDefDB,
      folder_name::String;
      replace_old::Bool=false,
      ignore_conflict::Bool=true,
      quarantine::Bool=true,
    ) -> Nothing

Move all the entries from the folder `folder_name` in the database of
definite lattices reached by `db` to the `main` folder.

See `promote_entry_to_main_folder!` for further documentation.
"""
function promote_folder_to_main_folder!(
  db::ZZLatDefDB,
  folder_name::String;
  replace_old::Bool=false,
  ignore_conflict::Bool=true,
  quarantine::Bool=true,
)
  _path = joinpath(path(db), folder_name)
  !isdir(_path) && return nothing
  for entry_label in readdir(_path)
    promote_entry_to_main_folder!(db, folder_name, entry_label; replace_old, ignore_conflict, quarantine)
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
      V::Vector{Int},
      n::Int,
      s::QQFieldElem = QQ(1),
    ) -> ZZLat

Return the integer lattice ``L(s)`` where ``L`` has rank ``n`` and the upper
triangular part of the Gram matrix of ``L`` is given by the entries in `V`.
"""
function lattice_from_data(
  V::Vector{Int},
  n::Int,
  s::QQFieldElem=QQ(1),
)
  @req length(V) == binomial(n+1, 2) "Wrong number of entries for a Gram matrix of rank $n"
  M = zero_matrix(QQ, n, n)
  k = 0
  for i in 1:n
    for j in i:n
      k += 1
      M[i, j] = M[j, i] = isone(s) ? QQ(V[k]) : s * V[k]
    end
  end
  return integer_lattice(; gram=M, cached=false, check=false)
end

@doc raw"""
    load_genus(
      db::ZZLatDefDB,
      G::ZZGenus,
    ) -> Vector{ZZLat}

If ``G`` can be loaded from the `main` folder of the database of
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
  lats = load_genus(db, label, s)
  genus(first(lats)) == G && return lats
  # Entries written by `save_genus!` before it rescaled to the reduced genus
  # are stored in the scaling of the genus they were saved for (same label):
  # read off that scaling s0 from the stored lattices, L_stored = L_red(s0).
  _, s0 = get_label_and_scaling_factor(genus(first(load_genus(db, label))))
  lats = load_genus(db, label, s // s0)
  genus(first(lats)) == G && return lats
  throw(GenusError(joinpath(path(db), "main", label), "No rescaling of the stored lattices lies in $G"))
end

@doc raw"""
    load_genus(
      db::ZZLatDefDB,
      entry_label::String,
      s::QQFieldElem = QQ(1);
      folder_name::String="main",
    ) -> Vector{ZZLat}

Return the list of lattices ``L(s)`` where ``L`` ranges over all the lattices
stored at the entry `entry_label` from the folder `folder_name` in the database
of definite lattices reached by `db`.

The keyword argument `folder_name` refers to the `main` folder by
default, but one can also input any other folder which is not `quarantine`.

If the entry `entry_label` from the folder `folder_name` does not exist, is
empty or is corrupted, an error is thrown.
"""
function load_genus(
  db::ZZLatDefDB,
  entry_label::String,
  s::QQFieldElem = QQ(1);
  folder_name::String="main",
)
  n = _rank_from_label(entry_label)
  _path = joinpath(path(db), folder_name, entry_label)
  @req isfile(_path) "Entry $entry_label does not exist in the folder $folder_name"
  lats = ZZLat[load_lattice(l, n, s; path=_path) for l in _read_entry(_path)]
  isempty(lats) && throw(GenusError(_path, "Empty entry"))
  !check_mass_formula(lats) && throw(MassFormulaError(_path))
  return lats
end

# the (nonempty) lines of the entry stored at `file_path`
function _read_entry(file_path::String)
  content = try
              f = ZipFile.Reader(file_path)
              try
                read(first(f.files), String)
              finally
                close(f)
              end
            catch e
              e isa InterruptException && rethrow()
              throw(GenusError(file_path, "Cannot read the zip archive"))
            end
  return [String(l) for l in eachsplit(content, '\n') if !isempty(strip(l))]
end

# the first line of the entry stored at `file_path` (without reading the rest)
function _read_first_line(file_path::String)
  f = ZipFile.Reader(file_path)
  try
    return readline(first(f.files))
  finally
    close(f)
  end
end

@doc raw"""
    load_lattice(
      l::AbstractString,
      n::Int,
      s::QQFieldElem = QQ(1);
      path::String = "",
    ) -> ZZLat

Return the lattice ``L(s)`` described by the line `l` of an entry of the
database (see `data_from_lattice`), where ``L`` has rank ``n``. The order of
its automorphism group is set. `path` is only used in error messages.
"""
function load_lattice(
  l::AbstractString,
  n::Int,
  s::QQFieldElem = QQ(1);
  path::String = "",
)
  where_ = isempty(path) ? string(l) : "$path (line $(repr(l)))"
  if !isone(count(isequal('|'), l))
    throw(GenusError(where_, "Wrong line format"))
  end

  lat, ord = split(strip(l), "|")

  _, V = try
           Hecke._parse(Vector{Int}, IOBuffer(lat))
         catch
           throw(GenusError(where_, "Cannot parse half gram"))
         end

  L = try
        lattice_from_data(V, n, s)
      catch
        throw(GenusError(where_, "Wrong number of entries in the half gram"))
      end

  L.automorphism_group_order = try
                                 last(Hecke._parse(ZZRingElem, IOBuffer(ord)))
                               catch
                                 throw(GenusError(where_, "Cannot parse automorphism group order"))
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
nonnegative signature, return a string to be written as a line in a file in
which the genus of ``M`` is stored. The string describes 2 pieces of
information separated by a `|`:
- the first piece encodes the entries of the upper triangular part of the
  Gram matrix of ``M``,
- the second piece encodes the order of the automorphism group of ``L``
  (which is the same as for ``M``).
"""
function data_from_lattice(L::ZZLat, s::QQFieldElem = QQ(1))
  M = gram_matrix(L)
  io = IOBuffer()
  print(io, "[")
  for i in 1:nrows(M), j in i:ncols(M)
    (i, j) == (1, 1) || print(io, ",")
    z = isone(s) ? M[i, j] : M[i, j] // s
    @req is_integral(z) "The lattice is not the rescaling by $s of an integral lattice"
    print(io, numerator(z))
  end
  print(io, "]|", automorphism_group_order(L))
  return String(take!(io))
end

# write the lines into the zip archive `file_path` (one member named `label`)
function _write_entry(file_path::String, label::String, lines)
  w = ZipFile.Writer(file_path)
  try
    f = ZipFile.addfile(w, label; method=ZipFile.Deflate)
    write(f, join(lines, "\n"))
  finally
    close(w)
  end
  return nothing
end

# a new empty file in the temporary folder, named `<label>__<random>`
function _temporary_path(db::ZZLatDefDB, label::String)
  tmp_path, io = mktemp(joinpath(path(db), "temporary"); cleanup=false)
  close(io)
  new_path = joinpath(path(db), "temporary", label*"__"*basename(tmp_path))
  mv(tmp_path, new_path)
  return new_path
end

@doc raw"""
    save_genus!(
      db::ZZLatDefDB,
      lats::Vector{ZZLat},
      G::ZZGenus = genus(first(lats));
      safe::Bool=true,
      save_duplicates::Bool=false,
      verbose::Bool=false,
    ) -> Bool

Update the `main` folder of the database of definite lattice reached by
`db` by storing the content of `lats`, rescaled to the reduced genus of
``G``, in a new entry labeled by the reduced genus of ``G``. Other ways to
call this function:
- `setindex!(db, lats, G)`
- `db[G] = lats`

If `safe` is set to `true`, the function checks whether the data is safe to be
saved (highly recommended). It essentially checks:
  - whether all the lattices in `lats` are of the genus `G`,
  - whether the mass formula for `G` is satisfied by the lattices in `lats`.
If the data is not safe, it is stored in the `quarantine` folder.

If the genus ``G`` is already available in `main`, then:
  - either the corresponding entry in `main` is corrupted. In this case, it
    is moved to `quarantine` and replaced by the new entry.
  - or `save_duplicates` is set to `true` and the new entry is saved in the
    folder `duplicate`.
  - or the new entry is ignored.

Return whether the genus is available in the `main` folder afterwards, i.e.
`false` if the data went to `quarantine`.

If `verbose` is set to `true`, the user allows the function to tell which
decision is taken for storing the new entry.
"""
function save_genus!(
  db::ZZLatDefDB,
  lats::Vector{ZZLat},
  G::ZZGenus = genus(first(lats));
  safe::Bool=true,
  save_duplicates::Bool=false,
  verbose::Bool=false,
)
  @req !isempty(lats) "No lattices to save"
  label, s = get_label_and_scaling_factor(G)
  verbose && println("Store new genus: ", label)
  add_new_key = true
  main_path = joinpath(path(db), "main", label)
  _path = main_path

  if safe
    if !all(isequal(G)∘genus, lats) || !check_mass_formula(lats)
      add_new_key = false
      verbose && println("!! Data corrupted !! The lattices are not all in the genus or do not satisfy the mass formula...")
      quarantine_label = _free_label(db, "quarantine", label, "qua")
      _path = joinpath(path(db), "quarantine", quarantine_label)
      verbose && println(">>>> New data stored in quarantine folder: ", quarantine_label)
    end
  end

  if add_new_key && isfile(main_path) && !save_duplicates && !is_corrupted_entry(db, label)
    verbose && println(">>>> The genus is already in the database, new data is ignored")
    add_key!(db, label)
    return true
  end

  # Write to the temporary folder first, check, then move into place: the
  # entry appears at its final place in one step (mv), so that concurrent
  # readers never see half a file.
  tmp_path = _temporary_path(db, label)
  tmp_label = basename(tmp_path)
  # the lattices are stored as lattices of the reduced genus: L = M(s)
  _write_entry(tmp_path, label, String[data_from_lattice(L, s) for L in lats])

  if add_new_key && is_corrupted_entry(db, tmp_label; folder_name="temporary")
    add_new_key = false
    verbose && println("!! Data corrupted !! Something went wrong during saving...")
    quarantine_label = _free_label(db, "quarantine", label, "qua")
    _path = joinpath(path(db), "quarantine", quarantine_label)
    verbose && println(">>>> New data stored in quarantine folder: ", quarantine_label)
  end

  if add_new_key && isfile(main_path)
    verbose && println("!! Duplicate warning !! The genus seems to be already in the database...")
    if is_corrupted_entry(db, label)
      verbose && println(">> Current version of genus is corrupted...")
      quarantine_label = move_to_quarantine!(db, "main", label)
      verbose && println(">>>> Old data moved to quarantine folder: ", quarantine_label)
      verbose && println(">>>> New data stored in main folder: ", label)
    else
      # (only reached with save_duplicates = true)
      verbose && println(">> Current version of genus is not corrupted...")
      duplicate_label = _free_label(db, "duplicate", label, "dup")
      _path = joinpath(path(db), "duplicate", duplicate_label)
      verbose && println(">>>> New data stored in duplicate folder: ", duplicate_label)
      mv(tmp_path, _path)
      add_key!(db, label)
      return true
    end
  elseif add_new_key
    verbose && println(">>>> New data stored in main folder: ", label)
  end

  mv(tmp_path, _path; force = (_path == main_path))
  add_new_key && add_key!(db, label)
  return add_new_key
end

###############################################################################
#
#  Normalization of entries stored in the wrong scaling
#
###############################################################################

@doc raw"""
    normalize_entry!(db::ZZLatDefDB, label::String) -> Bool

Rewrite the entry `label` of the `main` folder of `db` in the scaling of the
reduced genus, if it is stored in another one (as `save_genus!` did before
version 0.1.0 when it was given a genus that is not reduced, e.g. a negative
definite one). Return whether the entry was rewritten.
"""
function normalize_entry!(db::ZZLatDefDB, label::String)
  file_path = joinpath(path(db), "main", label)
  @req isfile(file_path) "Entry $label does not exist in the main folder"
  n = _rank_from_label(label)
  _label, s0 = get_label_and_scaling_factor(genus(load_lattice(_read_first_line(file_path), n; path=file_path)))
  @req _label == label "The entry $label stores a genus with another label ($_label)"
  isone(s0) && return false
  lines = _read_entry(file_path)
  # The stored lattices are L_red(s0): divide the Gram entries by s0. This is
  # done on the text, the orders of the automorphism groups do not change.
  p, q = Int(numerator(s0)), Int(denominator(s0))   # entries v -> v*q/p
  new_lines = String[]
  io = IOBuffer()
  for l in lines
    lat, ord = split(strip(l), "|")
    _, V = try
             Hecke._parse(Vector{Int}, IOBuffer(lat))
           catch
             throw(GenusError(file_path, "Cannot parse half gram"))
           end
    print(io, "[")
    for (k, v) in enumerate(V)
      w, r = divrem(v * q, p)
      iszero(r) || throw(GenusError(file_path, "Entries not divisible by the scaling factor $s0"))
      k > 1 && print(io, ",")
      print(io, w)
    end
    print(io, "]|", strip(ord))
    push!(new_lines, String(take!(io)))
  end
  # check the result before it replaces the entry
  _label2, s1 = get_label_and_scaling_factor(genus(load_lattice(first(new_lines), n; path=file_path)))
  @req _label2 == label && isone(s1) "Rewriting $label failed, the entry is left as it is"
  tmp_path = _temporary_path(db, label)
  _write_entry(tmp_path, label, new_lines)
  @req _read_entry(tmp_path) == new_lines "Rewriting $label failed, the entry is left as it is"
  mv(tmp_path, file_path; force=true)
  return true
end

@doc raw"""
    normalize_database!(db::ZZLatDefDB; verbose::Bool=false) -> Int

Apply `normalize_entry!` to all entries of the `main` folder of `db`, and
return the number of rewritten entries.
"""
function normalize_database!(db::ZZLatDefDB; verbose::Bool=false)
  update_keys!(db)
  n = 0
  for label in sort!(collect(index(db)))
    if normalize_entry!(db, label)
      n += 1
      verbose && println("rewritten: ", label)
    end
  end
  return n
end

###############################################################################
#
#  Printing
#
###############################################################################

_quantity(n::Int, singular::String, plural::String) = "$n " * (isone(n) ? singular : plural)

function Base.show(io::IO, ::MIME"text/plain", s::DBKeys{T}) where T
  io = AA.pretty(io)
  println(io, "DBKeys{$T}")
  print(io, AA.Indent(), "indexing a collection with ", _quantity(length(index(s)), "entry", "entries"))
  print(io, AA.Dedent())
end

function Base.show(io::IO, s::DBKeys{T}) where T
  print(io, "DBKeys{$T}")
end

function Base.show(io::IO, ::MIME"text/plain", db::ZZLatDefDB)
  io = AA.pretty(io)
  println(io, "Database of definite lattices")
  print(io, AA.Indent())
  print(io, "with ", _quantity(length(index(db)), "(reduced) genus", "(reduced) genera"))
  print(io, AA.Dedent())
end

function Base.show(io::IO, db::ZZLatDefDB)
  print(io, "Database of definite lattices")
end

###############################################################################
#
#  Migration from the uncompressed format
#
###############################################################################

# Compress the entries of the old (plain text) format found in `source` into
# `dist`; entries already present in `dist` are skipped.
function compress_old_data(
  db::ZZLatDefDB,
  source::String, # path to source
  dist::String = joinpath(path(db), "main") # path to dist
)
  tg = Set(readdir(dist))
  for label in readdir(source)
    genus_path = joinpath(dist, label)
    if label in tg
      println("skipped $label")
      continue
    end
    println("writing: $label")
    _write_entry(genus_path, label, readlines(joinpath(source, label)))
  end
  return nothing
end
