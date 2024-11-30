import Base: show

struct Mesh
    # basic
    grid::Matrix{Float64}
    edge2grid::Matrix{Int}
    face2edge::Matrix{Int}
    elem2face::Matrix{Int}
    # number
    Ndims::Int
    Ngrid::Int
    Nedge::Int
    Nface::Int
    Nelem::Int
    # auxillary index with fixed length
    face2grid::Matrix{Int}
    elem2edge::Matrix{Int}
    elem2grid::Matrix{Int}
    # auxillary index with various length
    grid2edge::Vector{Vector{Int}}
    grid2face::Vector{Vector{Int}}
    grid2elem::Vector{Vector{Int}}
    edge2face::Vector{Vector{Int}}
    edge2elem::Vector{Vector{Int}}
    face2elem::Vector{Vector{Int}}
    # normal Vector
    normx::Matrix{Float64}
    normy::Matrix{Float64}
    normz::Matrix{Float64}
end

"""
"""
function Mesh(grid::Matrix{Float64}, edge2grid::Matrix{Int}, face2edge::Matrix{Int}, elem2face::Matrix{Int})
    Ndims = size(grid, 1)
    Ngrid = size(grid, 2)
    Nedge = size(edge2grid, 2)
    Nface = size(face2edge, 2)
    Nelem = size(elem2face, 2)

    # allocate
    face2grid = zeros(Int, 3, Nface) #
    elem2edge = zeros(Int, 6, Nelem) #
    elem2grid = zeros(Int, 4, Nelem)

    grid2edge = Vector{Vector{Int}}(undef, Ngrid) #
    grid2face = Vector{Vector{Int}}(undef, Ngrid) #
    grid2elem = Vector{Vector{Int}}(undef, Ngrid) #
    for i = 1:Ngrid
        grid2edge[i] = Int[]
        grid2face[i] = Int[]
        grid2elem[i] = Int[]
    end
    edge2face = Vector{Vector{Int}}(undef, Nedge) #
    edge2elem = Vector{Vector{Int}}(undef, Nedge) #
    for i = 1:Nedge
        edge2face[i] = Int[]
        edge2elem[i] = Int[]
    end
    face2elem = Vector{Vector{Int}}(undef, Nface) #
    for i = 1:Nface
        face2elem[i] = Int[]
    end

    normx = zeros(size(elem2face, 1), Nelem)
    normy = zeros(size(elem2face, 1), Nelem)
    normz = zeros(size(elem2face, 1), Nelem)

    # build index
    for iedge = axes(edge2grid, 2)
        for ip = axes(edge2grid, 1)
            local ig = edge2grid[ip, iedge]
            if !in(iedge, grid2edge[ig])
                push!(grid2edge[ig], iedge)
            end
        end
    end

    for iface = axes(face2edge, 2)
        for iel = axes(face2edge, 1)
            local iedge = abs(face2edge[iel, iface])
            if !in(iface, edge2face[iedge])
                push!(edge2face[iedge], iface)
            end
            for igl = axes(edge2grid, 1)
                local igrid = edge2grid[igl, iedge]
                if !in(igrid, face2grid[:, iface])
                    jgl = findfirst(iszero, face2grid[:, iface])
                    face2grid[jgl, iface] = igrid
                end
                if !in(iface, grid2face[igrid])
                    push!(grid2face[igrid], iface)
                end
            end
        end
    end

    for ielem = axes(elem2face, 2)
        for ifl = axes(elem2face, 1)
            local iface = abs(elem2face[ifl, ielem])
            if !in(ielem, face2elem[iface])
                push!(face2elem[iface], ielem)
            end
            for iel = axes(face2edge, 1)
                local iedge = abs(face2edge[iel, iface])
                if !in(ielem, edge2elem[iedge])
                    push!(edge2elem[iedge], ielem)
                end
                if !in(iedge, elem2edge[:, ielem])
                    jel = findfirst(iszero, elem2edge[:, ielem])
                    elem2edge[jel, ielem] = iedge
                end
                for igl = axes(edge2grid, 1)
                    local igrid = edge2grid[igl, iedge]
                    if !in(ielem, grid2elem[igrid])
                        push!(grid2elem[igrid], ielem)
                    end
                    if !in(igrid, elem2grid[:, ielem])
                        jgl = findfirst(iszero, elem2grid[:, ielem])
                        elem2grid[jgl, ielem] = igrid
                    end
                end
            end
        end
    end

    for ielem = axes(elem2face, 2)
        for jface = axes(elem2face, 1)
            local iface = abs(elem2face[jface, ielem])
            fgrids = zeros(Int, size(elem2grid, 1))
            for jgrid = axes(face2grid, 1)
                fgrids[jgrid] = abs(face2grid[jgrid, iface])
            end
            for jgrid = axes(elem2grid, 1)
                local igrid = abs(elem2grid[jgrid, ielem])
                if !in(igrid, fgrids)
                    fgrids[end] = igrid
                end
            end
            r2 = grid[:, fgrids[2]] - grid[:, fgrids[1]]
            r3 = grid[:, fgrids[3]] - grid[:, fgrids[1]]
            r4 = grid[:, fgrids[4]] - grid[:, fgrids[1]]
            nm = normalize(cross(r2, r3))
            if dot(nm, r4) > 0.0
                nm = -nm
            end
            normx[jface, ielem] = nm[1]
            normy[jface, ielem] = nm[2]
            normz[jface, ielem] = nm[3]
        end
    end

    return Mesh(grid, edge2grid, face2edge, elem2face,
        Ndims, Ngrid, Nedge, Nface, Nelem,
        face2grid, elem2edge, elem2grid, grid2edge,
        grid2face, grid2elem, edge2face, edge2elem, face2elem,
        normx, normy, normz)
end

function show(io::IO, ::MIME"text/plain", m::Mesh)
    println(io, "Grid:")
    println(io, m.grid)
    println(io, "Edge:")
    println(io, m.edge2grid)
    println(io, "Face:")
    println(io, m.face2edge)
    println(io, "Element:")
    println(io, m.elem2face)
end