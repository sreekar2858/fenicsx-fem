"""
DOLFINx Modular Structural Solver - Cantilever beam test

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
# --------------------------------------
# 1. Imports & Config
# --------------------------------------
import yaml
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
from dolfinx import mesh, fem, io
import ufl
from slepc4py import SLEPc
from dolfinx.fem.petsc import assemble_matrix

# Load configuration
def load_config():
    with open("config.yaml") as f:
        return yaml.safe_load(f)
        
config = load_config()
material = config["material"]
geometry = config["geometry"]
mesh_params = config["mesh"]
solver_params = config["solver"]

# --------------------------------------
# 2. Mesh Generation
# --------------------------------------
def scale_stl(input_path, output_path, scale):
    """Scale an ASCII STL file by a given factor and write to output_path."""
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

def create_mesh():
    import gmsh
    from dolfinx.io import gmshio
    import os
    scale = geometry.get("scale", 1.0)
    # Option to use built-in BoxMesh for debugging mesh/solver pipeline
    if geometry.get("use_boxmesh", False):
        # Use BoxMesh like legacy FEniCS code
        L = geometry.get("L", 20.0)
        B = geometry.get("B", 0.5)
        H = geometry.get("H", 1.0)
        Nx = mesh_params.get("Nx", 200)
        Ny = mesh_params.get("Ny", int(B/L*Nx)+1)
        Nz = mesh_params.get("Nz", int(H/L*Nx)+1)
        domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0.,0.,0.]), np.array([L,B,H])], [Nx, Ny, Nz], cell_type=mesh.CellType.hexahedron)
        print("[DEBUG] Using built-in BoxMesh for domain.")
    else:
        gmsh.initialize([], False, False)  # Suppress gmsh terminal output
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cantilever")
        ext = os.path.splitext(geometry["input_file"])[1].lower()
        stl_file = geometry["input_file"]
        if ext == ".stl":
            if scale != 1.0:
                scaled_stl = stl_file.replace('.stl', f'_scaled_{scale}.stl')
                scale_stl(stl_file, scaled_stl, scale)
                stl_file = scaled_stl
            # Refine mesh by setting smaller mesh size
            gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
            gmsh.merge(stl_file)
            gmsh.model.mesh.classifySurfaces(40 * (3.14159 / 180), True, True, 180)
            gmsh.model.mesh.createGeometry()
            surfaces = gmsh.model.getEntities(dim=2)
            surface_tags = [s[1] for s in surfaces]
            sl = gmsh.model.geo.addSurfaceLoop(surface_tags)
            vol = gmsh.model.geo.addVolume([sl])
            gmsh.model.geo.synchronize()
            gmsh.model.addPhysicalGroup(3, [vol], tag=1)
        elif ext in [".step", ".stp", ".x_t", ".x_b"]:
            gmsh.model.occ.importShapes(geometry["input_file"])
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
        gmsh.model.mesh.setOrder(mesh_params["element_order"])
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
        from dolfinx.cpp.mesh import cell_volume
        vols = cell_volume(domain)
        if np.any(vols <= 0):
            print("[WARNING] Mesh has cells with non-positive volume!")
    except Exception:
        pass
    return domain

def create_boxmesh():
    # Match GMSH mesh: x=[-0.25,0.25], y=[0,20], z=[-0.5,0.5]
    x0, x1 = -0.25, 0.25
    y0, y1 = 0.0, 20.0
    z0, z1 = -0.5, 0.5
    Nx = mesh_params.get("Nx", 5)  # Adjust for similar resolution
    Ny = mesh_params.get("Ny", 100)
    Nz = mesh_params.get("Nz", 10)
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
    return domain

# Choose mesh type for main analysis
use_boxmesh_for_main = False  # Set to True to use BoxMesh for main, False for GMSH/STL
domain = create_boxmesh() if use_boxmesh_for_main else create_mesh()
# Export mesh to VTU for visualization (use VTKFile for .vtu extension)
with io.VTKFile(domain.comm, "mesh_export.vtu", "w") as vtk:
    vtk.write_mesh(domain)
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

fixed_cells = np.unique(fixed_facets)
bc_marker.x.array[fixed_cells] = 1

with io.VTKFile(domain.comm, "bc_marker.vtu", "w") as vtk:
    vtk.write_mesh(domain)
    vtk.write_function(bc_marker)

# --------------------------------------
# 3. Function Space & BCs
# --------------------------------------
V = fem.functionspace(domain, ("Lagrange", mesh_params["element_order"], (3,)))

def min_y_face(x):
    y_min = np.min(domain.geometry.x[:, 1])
    return np.isclose(x[1], y_min, atol=1e-6)

boundary_dofs = fem.locate_dofs_geometrical(V, min_y_face)
print(f"Number of fixed DOFs: {len(boundary_dofs)}")
u_bc = fem.Constant(domain, PETSc.ScalarType((0.0, 0.0, 0.0)))
bc = fem.dirichletbc(u_bc, boundary_dofs, V)

# --------------------------------------
# 4. Material Model
# --------------------------------------
E = float(material["E"])
nu = float(material["nu"])
rho = float(material["density"])
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
def epsilon(u):
    return ufl.sym(ufl.grad(u))
def sigma(u):
    return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(3) + 2.0 * mu * epsilon(u)

# --------------------------------------
# 5. Variational Problem
# --------------------------------------
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
m = rho * ufl.inner(u, v) * ufl.dx

# --------------------------------------
# 6. Assembly
# --------------------------------------
# Assemble stiffness matrix with BCs (mimics legacy FEniCS assemble_system)
# Note: DOLFINx does not use a dummy linear form, but this matches legacy intent
A = assemble_matrix(fem.form(a), bcs=[bc])
A.assemble()
M = assemble_matrix(fem.form(m), bcs=[])
M.assemble()

# --------------------------------------
# 7. Eigenvalue Solver
# --------------------------------------
num_modes = int(solver_params["num_modes"])
tol = float(solver_params["tol"])
max_it = int(solver_params["max_it"])
eigensolver = SLEPc.EPS().create(domain.comm)
eigensolver.setOperators(A, M)
eigensolver.setProblemType(SLEPc.EPS.ProblemType.GHEP)
eigensolver.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL)
eigensolver.setDimensions(num_modes, PETSc.DECIDE)
eigensolver.setTolerances(tol, max_it)
# Use shift-and-invert as in legacy FEniCS code
st = eigensolver.getST()
st.setType(SLEPc.ST.Type.SINVERT)
eigensolver.setTarget(0.0)  # Match legacy code shift
# eigensolver.setFromOptions()  # Only use if you want to override with command-line options
eigensolver.solve()
nconv = eigensolver.getConverged()
print(f"Converged {nconv} eigenpairs")
print(f"SLEPc reason: {eigensolver.getConvergedReason()}")

# Run modal analysis on BoxMesh for comparison (only if main mesh is not BoxMesh)
if not use_boxmesh_for_main:
    print("Runnning BoxMesh Case")
    box_domain = create_boxmesh()
    V_box = fem.functionspace(box_domain, ("Lagrange", mesh_params["element_order"], (3,)))
    # Apply BC on min-x face for BoxMesh (to match GMSH mesh orientation)
    def min_y_face_box(x):
        y_min = np.min(box_domain.geometry.x[:, 1])
        return np.isclose(x[1], y_min, atol=1e-6)
    boundary_dofs_box = fem.locate_dofs_geometrical(V_box, min_y_face_box)
    u_bc_box = fem.Constant(box_domain, PETSc.ScalarType((0.0, 0.0, 0.0)))
    bc_box = fem.dirichletbc(u_bc_box, boundary_dofs_box, V_box)
    E_box = float(material["E"])
    nu_box = float(material["nu"])
    rho_box = float(material["density"])
    mu_box = E_box / (2.0 * (1.0 + nu_box))
    lmbda_box = E_box * nu_box / ((1.0 + nu_box) * (1.0 - 2.0 * nu_box))
    def epsilon_box(u):
        return ufl.sym(ufl.grad(u))
    def sigma_box(u):
        return lmbda_box * ufl.tr(epsilon_box(u)) * ufl.Identity(3) + 2.0 * mu_box * epsilon_box(u)
    u_box = ufl.TrialFunction(V_box)
    v_box = ufl.TestFunction(V_box)
    a_box = ufl.inner(sigma_box(u_box), epsilon_box(v_box)) * ufl.dx
    m_box = rho_box * ufl.inner(u_box, v_box) * ufl.dx
    A_box = assemble_matrix(fem.form(a_box), bcs=[bc_box])
    A_box.assemble()
    M_box = assemble_matrix(fem.form(m_box), bcs=[])
    M_box.assemble()
    num_modes_box = int(solver_params["num_modes"])
    tol_box = float(solver_params["tol"])
    max_it_box = int(solver_params["max_it"])
    eigensolver_box = SLEPc.EPS().create(box_domain.comm)
    eigensolver_box.setOperators(A_box, M_box)
    eigensolver_box.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    eigensolver_box.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL)
    eigensolver_box.setDimensions(num_modes_box, PETSc.DECIDE)
    eigensolver_box.setTolerances(tol_box, max_it_box)
    st_box = eigensolver_box.getST()
    st_box.setType(SLEPc.ST.Type.SINVERT)
    eigensolver_box.setTarget(0.0)
    eigensolver_box.solve()
    nconv_box = eigensolver_box.getConverged()
    print(f"[BoxMesh] Converged {nconv_box} eigenpairs")
    for i in range(min(num_modes_box, nconv_box)):
        eigval = eigensolver_box.getEigenvalue(i)
        freq = np.sqrt(eigval.real) / (2 * np.pi) if eigval.real > 0 else float('nan')
        print(f"[BoxMesh] Mode {i+1}: λ = {eigval.real:.6e}, f = {freq:.6f} Hz")

# --------------------------------------
# 8. Output Results
# --------------------------------------
for i in range(min(num_modes, nconv)):
    eigval = eigensolver.getEigenvalue(i)
    freq = np.sqrt(eigval.real) / (2 * np.pi) if eigval.real > 0 else float('nan')
    print(f"Mode {i+1}: λ = {eigval.real:.6e}, f = {freq:.6f} Hz")
    eigvec = PETSc.Vec().createWithArray(np.zeros(V.dofmap.index_map.size_local * V.dofmap.index_map_bs, dtype=PETSc.ScalarType), comm=domain.comm)
    eigensolver.getEigenvector(i, eigvec)
    mode = fem.Function(V, name=f"Mode_{i+1}")
    mode.x.array[:] = eigvec.getArray().real
    mode.x.scatter_forward()
    with io.XDMFFile(domain.comm, f"mode_{i+1}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(mode)
print("\nModal analysis complete. Results written to mode_*.xdmf files.")



