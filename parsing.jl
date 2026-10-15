# Backwards compatibility: scripts that do
#
#     include("path/to/DefLatDB/parsing.jl")
#     db = definite_lattice_database("path/to/DefLatDB/data")
#
# keep working. They need Hecke and ZipFile in their environment (as before).
# New code should use the package instead:
#
#     using DefLatDB
#     db = definite_lattice_database()
if !isdefined(@__MODULE__, :DefLatDB)
  include(joinpath(@__DIR__, "src", "DefLatDB.jl"))
end
using .DefLatDB
# `path` is not exported by the package, but older scripts use it; `index` is
# a method of `Hecke.index`, available wherever Hecke (or Oscar) is loaded
if !isdefined(@__MODULE__, :path)
  using .DefLatDB: path
end
