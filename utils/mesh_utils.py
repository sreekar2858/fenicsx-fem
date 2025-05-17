"""
DOLFINx Modular Structural Solver - Mesh utilities

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
def scale_stl(input_path, output_path, scale):
    with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.strip().startswith('vertex'):
                parts = line.split()
                try:
                    nums = [float(x) for x in parts[1:]]
                    scaled = [str(float(x) * scale) for x in parts[1:]]
                    fout.write(f"{parts[0]} {' '.join(scaled)}\n")
                except Exception:
                    fout.write(line)
            elif line.strip().startswith('facet normal'):
                parts = line.split()
                try:
                    nums = [float(x) for x in parts[-3:]]
                    scaled = [str(float(x) * scale) for x in parts[-3:]]
                    fout.write(f"{' '.join(parts[:-3])} {' '.join(scaled)}\n")
                except Exception:
                    fout.write(line)
            else:
                fout.write(line)

# Create a mesh using GMSH or BoxMesh
def create_mesh(geometry_params=None, mesh_params=None, solver_params=None):
    import gmsh
    from dolfinx.io import gmshio
    import os
    import numpy as np
    from mpi4py import MPI
    if geometry_params.get("use_boxmesh", False):
        from dolfinx import mesh
        # Use BoxMesh like legacy FEniCS code
        x0 = geometry_params.get("dimensions").get("x")[0]
        y0 = geometry_params.get("dimensions").get("y")[0]
        z0 = geometry_params.get("dimensions").get("z")[0]
        x1 = geometry_params.get("dimensions").get("x")[1]
        y1 = geometry_params.get("dimensions").get("y")[1]
        z1 = geometry_params.get("dimensions").get("z")[1]
        Nx = mesh_params.get("divisions").get("Nx", None)
        Ny = mesh_params.get("divisions").get("Ny", None)
        Nz = mesh_params.get("divisions").get("Nz", None)
        domain = mesh.create_box(
                MPI.COMM_WORLD,
                [np.array([x0, y0, z0]), np.array([x1, y1, z1])],
                [Nx, Ny, Nz],
                cell_type=mesh.CellType.hexahedron
            )
        print("[BoxMesh] Using built-in BoxMesh for domain.")
        num_cells = domain.topology.index_map(domain.topology.dim).size_local
        num_vertices = domain.topology.index_map(0).size_local
        print(f"[BoxMesh] Mesh: {num_cells} cells, {num_vertices} vertices")
        print(f"[BoxMesh] Mesh geometric dimension: {domain.geometry.dim}")
        print(f"[BoxMesh] Mesh coordinate shape: {domain.geometry.x.shape}")
        print(f"[BoxMesh] Mesh coordinate min: {np.min(domain.geometry.x, axis=0)}")
        print(f"[BoxMesh] Mesh coordinate max: {np.max(domain.geometry.x, axis=0)}")        
    else:
        gmsh.initialize([], False, False)  # Suppress gmsh terminal output
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(solver_params.get("name", "GMSH"))
        scale = float(geometry_params.get("scale", 1.0))
        ext = os.path.splitext(geometry_params.get("input_file"))[1].lower()
        stl_file = geometry_params.get("input_file")
        if ext == ".stl":
            if scale != 1.0:
                scaled_stl = stl_file.replace('.stl', f'_scaled_{scale}.stl')
                scale_stl(stl_file, scaled_stl, scale)
                stl_file = scaled_stl
            # Refine mesh by setting smaller mesh size
            if mesh_params.get("min_size", None):
                gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_params.get("min_size"))
            if mesh_params.get("max_size", None):
                gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_params.get("max_size"))
            gmsh.merge(stl_file)
            gmsh.model.mesh.classifySurfaces(40 * (3.14159 / 180), True, True, 180)
            gmsh.model.mesh.createGeometry()
            surfaces = gmsh.model.getEntities(dim=2)
            surface_tags = [s[1] for s in surfaces]
            sl = gmsh.model.geo.addSurfaceLoop(surface_tags)
            vol = gmsh.model.geo.addVolume([sl])
            gmsh.model.geo.synchronize()
            gmsh.model.addPhysicalGroup(3, [vol], tag=1)
        # TODO: Need to add support for other file formats UNTESTED
        elif ext in [".step", ".stp", ".x_t", ".x_b"]:
            gmsh.model.occ.importShapes(geometry_params["input_file"])
            if scale != 1.0:
                gmsh.model.occ.synchronize()
                all_entities = gmsh.model.occ.getEntities(3) or gmsh.model.occ.getEntities(2)
                for dim, tag in all_entities:
                    gmsh.model.occ.scale([(dim, tag)], 0, 0, 0, scale, scale, scale)
                gmsh.model.occ.synchronize()
            volumes = gmsh.model.getEntities(dim=3)
            volume_tags = [v[1] for v in volumes]
            if volume_tags:
                gmsh.model.addPhysicalGroup(3, volume_tags, tag=1)
            else:
                raise RuntimeError("No volumes found in CAD file.")
        else:
            raise RuntimeError(f"Unsupported file extension: {ext}")
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(mesh_params.get("element_order", 1))
        from dolfinx.io import gmshio
        domain, _, _ = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
        gmsh.finalize()
    # Print mesh stats
    num_cells = domain.topology.index_map(domain.topology.dim).size_local
    num_vertices = domain.topology.index_map(0).size_local
    print(f"Mesh: {num_cells} cells, {num_vertices} vertices")
    print(f"Mesh geometric dimension: {domain.geometry.dim}")
    print(f"Mesh coordinate shape: {domain.geometry.x.shape}")
    print(f"Mesh coordinate min: {np.min(domain.geometry.x, axis=0)}")
    print(f"Mesh coordinate max: {np.max(domain.geometry.x, axis=0)}")
    # Mesh validity checks
    if domain.topology.dim != 3:
        raise RuntimeError(f"Mesh is not 3D! dim={domain.topology.dim}")
    if num_cells == 0:
        raise RuntimeError("Mesh has zero cells!")
    # Check for degenerate cells (zero or negative volume)
    try:
        # WARNING: The import doesn't work
        from dolfinx.cpp.mesh import cell_volume
        vols = cell_volume(domain)
        if np.any(vols <= 0):
            print("[WARNING] Mesh has cells with non-positive volume!")
    except Exception:
        pass
    return domain

# In utils/mesh_utils.py
def mark_boundary_facets(domain):
    from dolfinx import mesh, fem
    import numpy as np
    # WARNING: This function is hardcoded to mark the fixed face (x = min face)
    # Mark boundary condition face for ParaView visualization
    marker_space = fem.functionspace(domain, ("DG", 0))
    bc_marker = fem.Function(marker_space)
    bc_marker.x.array[:] = 0

    # Find facets on the fixed face (x = min face, to match legacy FEniCS)
    fdim = domain.topology.dim - 1
    domain.topology.create_connectivity(fdim, domain.topology.dim)
    exterior_facets = mesh.exterior_facet_indices(domain.topology)
    fixed_facets = []
    x_min = np.min(domain.geometry.x[:, 0])
    for f in exterior_facets:
        verts = domain.topology.connectivity(fdim, 0).links(f)
        coords = domain.geometry.x[verts]
        if np.all(np.isclose(coords[:, 0], x_min, atol=1e-6)):
            cells = domain.topology.connectivity(fdim, domain.topology.dim).links(f)
            fixed_facets.extend(cells)

    marked_cells = np.unique(fixed_facets)
    bc_marker.x.array[marked_cells] = 1
    
    return bc_marker, marked_cells