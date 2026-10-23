# DefLatDB

A database of genera of definite integer lattices. Every entry is a **complete
list of representatives of the isometry classes** of a genus, together with
the orders of their automorphism groups. The data was computed with
[Hecke](https://github.com/thofma/Hecke.jl) and
[Oscar](https://github.com/oscar-system/Oscar.jl) (Kneser's neighbour method,
isometries of prime order, root sublattices). It is checked against the mass
formula.

This repository is a Julia package to read and extend the database. Its only
dependencies are Hecke and ZipFile (for the compressed entries).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/StevellM/DefLatDB.git")
```

The data (`data/main`, about 0.8 GB) is part of the repository.

## Usage

```julia
using Hecke, DefLatDB

db = definite_lattice_database()       # the data shipped with the package
# db = definite_lattice_database("/path/to/DefLatDB/data")   # another copy

G = genus(root_lattice(:E, 8))
haskey(db, G)          # is the genus in the database?
R = db[G]              # the representatives (Vector{ZZLat}) ...
automorphism_group_order(R[1])   # ... with the orders of their automorphism groups

# Any rescaling of a stored genus is found as well, e.g. negative definite ones:
db[rescale(G, -1)]

# Store a genus (after checking that `lats` lies in G and satisfies the mass formula):
save_genus!(db, lats, G)          # or db[G] = lats
```

`db[G]` returns lattices in the genus `G` itself. One entry serves all
rescalings `G(s)` of a genus.

## How the data is stored

```
data/
  main/         the database: one file per (reduced) genus
  duplicate/    second copies of genera already in main (save_genus!(...; save_duplicates=true))
  quarantine/   data that failed a check
  temporary/    files being written
```

Only `data/main` is part of the repository; the other folders are created on
demand.

- **Reduced genus.** Each genus is stored once, as its *reduced genus*: the
  rescaling that is integral, primitive (scale 1) and of nonnegative signature.
  `get_label_and_scaling_factor(G)` returns its label and the factor `s` with
  `G = G_red(s)`.
- **Label = file name.** The label is `r.p.d.l.t` followed by the genus symbol of the
  reduced genus: rank, positive index of inertia, determinant, level, type
  (`I` odd, `II` even). For example, `9.9.28.28.II{1}^{-8}_{II}{4}^{-1}_{3}{7}^{-1}`
  is the genus `II_(9, 0) 4^-1_3 7^-1`.
- **Entry.** Each file is a zip archive with one member. Its lines describe one
  lattice each: `[g11,g12,...,g1n,g22,...,gnn]|o`, i.e. the upper triangle of a
  Gram matrix (row by row) and the order of the automorphism group.

## Checks and maintenance

```julia
is_valid_entry(db, label)         # same genus, pairwise non-isometric, mass formula (slow)
test_validity_database(db)        # is_valid_entry for all entries; returns the number of invalid ones
cleanup_database!(db)             # move entries that fail the quick checks to quarantine
promote_folder_to_main_folder!(db, "duplicate")
empty_folder!(db, "quarantine")
update_keys!(db)                  # re-read the index (entries stored by other processes)
normalize_database!(db)           # rewrite entries that are not in the scaling of the reduced genus
```

`load_genus` raises a `GenusError` or `MassFormulaError` for a corrupted entry.

Entries are written to `temporary/` and moved into place in one step, so
several processes may share one database directory.

## Old scripts

Scripts that load the database with

```julia
include("path/to/DefLatDB/parsing.jl")
db = definite_lattice_database("path/to/DefLatDB/data")
```

keep working. Their environment needs Hecke (or Oscar) and ZipFile, as before.

## Tests

```julia
using Pkg; Pkg.test("DefLatDB")
```

## Citing us

If you use the database in a substantial way for your research, please consider
citing us in the following way:

```bib
@misc{DefLatDB,
  author = {Brandhorst, Simon and Muller, Stevell},
  title = {DefLatDB: A database of genera of definite integer lattices, {V}ersion 0.1.0},
  year = {2026},
  url = {https://github.com/StevellM/DefLatDB},
}
```
