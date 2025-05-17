"""
DOLFINx Modular Structural Solver - Core solver logic

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import assemble_matrix
from petsc4py import PETSc
from slepc4py import SLEPc
import ufl

import os
import numpy as np
from utils.file_utils import load_config
from utils.mesh_utils import (create_mesh, mark_boundary_facets)

def run_modal_analysis(config_path):
    # --- Load configuration ---
    config = load_config(config_path)

    material_params = config["material"]
    geometry_params = config["geometry"]
    mesh_params = config["mesh"]
    solver_params = config["solver"]
    export_params = config["exports"]

    # --- Create mesh ---
    domain = create_mesh(geometry_params, mesh_params, solver_params)
    if export_params.get("mesh", None):
        with io.VTKFile(domain.comm, os.path.join("output", export_params["mesh"]["filename"]), "w") as vtk:
            vtk.write_mesh(domain)

    # --- Dirichlet boundary condition ---
    def min_y_face(x):
        y_min = np.min(domain.geometry.x[:, 1])
        return np.isclose(x[1], y_min, atol=1e-6)

    V = fem.functionspace(domain, ("Lagrange", 1 if geometry_params.get("use_boxmesh") else mesh_params.get("element_order", 1), (3,)))
    boundary_dofs = fem.locate_dofs_geometrical(V, min_y_face)
    print(f"Number of fixed DOFs: {len(boundary_dofs)}")
    u_bc = fem.Constant(domain, PETSc.ScalarType((0.0, 0.0, 0.0)))
    bc = fem.dirichletbc(u_bc, boundary_dofs, V)

    # Mark the boundary
    bc_marker, fixed_cells = mark_boundary_facets(domain)
    if export_params.get("bc_marker", None):
        with io.VTKFile(domain.comm, os.path.join("output", export_params["bc_marker"]["filename"]), "w") as vtk:
            vtk.write_mesh(domain)
            vtk.write_function(bc_marker)

    # --- Material properties ---
    E = float(material_params["E"])
    nu = float(material_params["nu"])
    rho = float(material_params["density"])
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    def epsilon(u):
        return ufl.sym(ufl.grad(u))
    def sigma(u):
        return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(3) + 2.0 * mu * epsilon(u)

    # --- Variational problem ---
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    m = rho * ufl.inner(u, v) * ufl.dx

    # --- Assemble system ---
    A = assemble_matrix(fem.form(a), bcs=[bc])
    A.assemble()
    M = assemble_matrix(fem.form(m), bcs=[])
    M.assemble()

    # --- Solve system ---
    num_modes = int(solver_params.get("num_modes", 6))
    tol = float(solver_params.get("tol", 1e-8))
    max_it = int(solver_params.get("max_it", 100))
    eigensolver = SLEPc.EPS().create(domain.comm)
    eigensolver.setOperators(A, M)
    eigensolver.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    eigensolver.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL)
    eigensolver.setDimensions(num_modes, PETSc.DECIDE)
    eigensolver.setTolerances(tol, max_it)
    st = eigensolver.getST()
    st.setType(SLEPc.ST.Type.SINVERT)
    eigensolver.setTarget(0.0)
    eigensolver.solve()
    nconv = eigensolver.getConverged()
    print(f"Converged {nconv} eigenpairs")
    print(f"SLEPc reason: {eigensolver.getConvergedReason()}")

    # --- Extract eigenvalues and eigenvectors ---
    for i in range(min(num_modes, nconv)):
        eigval = eigensolver.getEigenvalue(i)
        freq = np.sqrt(eigval.real) / (2 * np.pi) if eigval.real > 0 else float('nan')
        print(f"Mode {i+1}: λ = {eigval.real:.6e}, f = {freq:.6f} Hz")
        eigvec = PETSc.Vec().createWithArray(
            np.zeros(V.dofmap.index_map.size_local * V.dofmap.index_map_bs, dtype=PETSc.ScalarType),
            comm=domain.comm)
        eigensolver.getEigenvector(i, eigvec)
        mode = fem.Function(V, name=f"Mode_{i+1}")
        mode.x.array[:] = eigvec.getArray().real
        mode.x.scatter_forward()
        if export_params.get("modes", None):
            with io.XDMFFile(domain.comm, os.path.join("output", f"{export_params['modes']['filename'].split('.')[0]}_{i+1}.xdmf"), "w") as xdmf:
                xdmf.write_mesh(domain)
                xdmf.write_function(mode)
