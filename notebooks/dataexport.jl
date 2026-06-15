# assumes directory is available
function exportnum(tuple)
    datadir = "/home/gert/Projects/FusionRings.jl/src/data/split_number_data/"

    path = joinpath(datadir,tuple[1]*".mrdi")
    Oscar.save(path,tuple[2])
end
