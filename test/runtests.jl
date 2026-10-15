using Test
using Hecke
using DefLatDB
import DefLatDB: ZipFile, index, path

# a genus of 4 classes: II_(9, 0) 4^-1_3 7^-1
const G = genus(direct_sum(root_lattice(:A, 6), root_lattice(:A, 3))[1])
const R = enumerate_definite_genus(G)
const LABEL = first(get_label_and_scaling_factor(G))

# overwrite (or create) an entry with the given lines
function write_entry(db, folder, label, lines)
  w = ZipFile.Writer(joinpath(path(db), folder, label))
  f = ZipFile.addfile(w, label; method=ZipFile.Deflate)
  write(f, join(lines, "\n"))
  close(w)
end

@testset "DefLatDB" begin
  @testset "create, save, load" begin
    db = definite_lattice_database(mktempdir())
    for folder in ("main", "buffer", "quarantine", "duplicate", "temporary")
      @test isdir(joinpath(path(db), folder))
    end
    @test length(db) == 0
    @test !haskey(db, G)

    @test save_genus!(db, R, G)
    @test haskey(db, G)
    @test G in db
    @test haskey(db, LABEL)
    @test isempty(readdir(joinpath(path(db), "temporary")))
    R2 = db[G]
    @test length(R2) == length(R)
    @test all(L -> genus(L) == G, R2)
    @test sum(L -> 1//automorphism_group_order(L), R2) == mass(G)
    @test length(db[LABEL]) == length(R)

    # the same entry serves all rescalings of the genus
    for c in (-1, 2, -3//2)
      Gc = rescale(G, c)
      @test haskey(db, Gc)
      @test all(L -> genus(L) == Gc, db[Gc])
    end

    # a fresh object sees the entry on disk
    db2 = definite_lattice_database(path(db))
    @test haskey(db2, G)
    @test is_valid_entry(db2, LABEL)
    @test test_validity_database(db2) == 0
  end

  @testset "a rescaled genus is stored as the reduced genus" begin
    # (save_genus! used to store the lattices in the scaling of the genus they
    # were saved for, under the label of the reduced genus)
    db = definite_lattice_database(mktempdir())
    Gn = rescale(G, -2)
    @test save_genus!(db, ZZLat[rescale(L, -2) for L in R], Gn)
    @test all(L -> genus(L) == G, load_genus(db, LABEL))   # what is on disk
    @test all(L -> genus(L) == Gn, db[Gn])
    @test all(L -> genus(L) == G, db[G])

    # entries stored the old way are still loaded correctly ...
    write_entry(db, "main", LABEL, [data_from_lattice(rescale(L, -2)) for L in R])
    @test all(L -> genus(L) == Gn, db[Gn])
    @test all(L -> genus(L) == G, db[G])
    # ... and can be rewritten in the reduced scaling
    @test normalize_database!(db) == 1
    @test all(L -> genus(L) == G, load_genus(db, LABEL))
    @test [automorphism_group_order(L) for L in load_genus(db, LABEL)] == [automorphism_group_order(L) for L in R]
    @test !normalize_entry!(db, LABEL)
    @test isempty(readdir(joinpath(path(db), "temporary")))
  end

  @testset "duplicates, quarantine, buffer" begin
    db = definite_lattice_database(mktempdir())
    @test save_genus!(db, R, G)
    # already there: ignored, or stored as a duplicate
    @test save_genus!(db, R, G)
    @test isempty(readdir(joinpath(path(db), "duplicate")))
    @test save_genus!(db, R, G; save_duplicates=true)
    dups = readdir(joinpath(path(db), "duplicate"))
    @test dups == [LABEL * "__dup.0"]
    promote_folder_to_main_folder!(db, "duplicate")   # conflict: the duplicate is dropped
    @test isempty(readdir(joinpath(path(db), "duplicate")))
    @test haskey(db, G)

    # incomplete data goes to quarantine, under distinct names
    db3 = definite_lattice_database(mktempdir())
    @test !save_genus!(db3, R[1:end-1], G)
    @test !save_genus!(db3, R[1:end-1], G)
    @test sort(readdir(joinpath(path(db3), "quarantine"))) == [LABEL * "__qua.0", LABEL * "__qua.1"]
    @test !haskey(db3, G)
    @test isempty(readdir(joinpath(path(db3), "temporary")))

    # a corrupted entry in main is detected and replaced...
    write_entry(db3, "main", LABEL, [data_from_lattice(L) for L in R[1:end-1]])
    update_keys!(db3)
    @test haskey(db3, G)
    @test is_corrupted_entry(db3, LABEL)
    @test save_genus!(db3, R, G)
    @test !is_corrupted_entry(db3, LABEL)
    @test length(readdir(joinpath(path(db3), "quarantine"))) == 3

    # ... or the new data goes to the buffer, from where it can be promoted
    write_entry(db3, "main", LABEL, [data_from_lattice(L) for L in R[1:end-1]])
    @test !save_genus!(db3, R, G; replace_corrupted_data=false)
    @test !haskey(db3, G)                    # the corrupted entry was removed from the index
    @test readdir(joinpath(path(db3), "buffer")) == [LABEL * "__buf.0"]
    promote_folder_to_main_folder!(db3, "buffer")
    @test haskey(db3, G)                     # (the key is the label, not the buffer label)
    @test !haskey(db3, LABEL * "__buf.0")
    @test is_valid_entry(db3, LABEL)

    # cleanup moves corrupted entries out of main
    write_entry(db3, "main", LABEL, [data_from_lattice(L) for L in R[1:end-1]])
    cleanup_database!(db3)
    @test !haskey(db3, G)
    @test !isfile(joinpath(path(db3), "main", LABEL))

    empty_folder!(db3, "quarantine")
    @test isempty(readdir(joinpath(path(db3), "quarantine")))
    @test_throws ArgumentError empty_folder!(db3, "main")
  end

  @testset "parsing" begin
    L = R[1]
    l = data_from_lattice(L)
    L2 = load_lattice(l, rank(L))
    @test gram_matrix(L2) == gram_matrix(L)
    @test automorphism_group_order(L2) == automorphism_group_order(L)
    @test genus(load_lattice(l, rank(L), QQ(-2))) == rescale(G, -2)
    @test_throws GenusError load_lattice("[2]|abc", 1)
    @test_throws GenusError load_lattice("[2,1]|2", 1)
    @test_throws GenusError load_lattice("[2]", 1)
    @test_throws ArgumentError data_from_lattice(L, QQ(2))   # L is not L'(2) for an integral L'

    db = definite_lattice_database(mktempdir())
    write_entry(db, "main", LABEL, [data_from_lattice(L) for L in R])
    # a trailing newline is fine
    write_entry(db, "main", LABEL, [[data_from_lattice(L) for L in R]; ""])
    @test length(load_genus(db, LABEL)) == length(R)
  end

  @testset "printing" begin
    db = definite_lattice_database(mktempdir())
    @test sprint(show, MIME"text/plain"(), db) == "Database of definite lattices\n  with 0 (reduced) genera"
    save_genus!(db, R, G)
    @test occursin("with 1 (reduced) genus", sprint(show, MIME"text/plain"(), db))
    @test occursin("1 entry", sprint(show, MIME"text/plain"(), keys(db)))
    @test sprint(show, db) == "Database of definite lattices"
  end

  @testset "the data shipped with the package" begin
    main = joinpath(DefLatDB.DEFAULT_DATA_PATH, "main")
    if isdir(main) && !isempty(readdir(main))
      db = definite_lattice_database()
      @test length(db) > 0
      for label in sort(collect(index(db)); by = l -> filesize(joinpath(main, l)))[1:20]
        lats = load_genus(db, label)
        Gl = genus(first(lats))
        @test first(get_label_and_scaling_factor(Gl)) == label
        @test haskey(db, Gl)
        @test all(L -> genus(L) == Gl, db[Gl])
      end
    end
  end
end
