# DOLFINx Modular Structural Solver
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![DOLFINx](https://img.shields.io/badge/dolfinx-0.1.0-blue.svg)](https://www.fenicsproject.org/download/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/dolfinx/badge/?version=latest)](https://dolfinx.readthedocs.io/en/latest/?badge=latest)
[![Code Coverage](https://codecov.io/gh/sreekar2858/dolfinx/branch/main/graph/badge.svg)](https://codecov.io/gh/sreekar2858/dolfinx)

<p align="center">
  <img src="https://fenicsproject.org/assets/logo/fenics_logo.svg" alt="Logo" width="100"/>
</p>

A modular, professional codebase for structural modal analysis using [DOLFINx](https://github.com/FEniCS/dolfinx), designed for clarity, flexibility, and reproducibility. This solver computes the natural frequencies and mode shapes of structures using the finite element method.

## Author

**Sreekar Reddy, Sajjala**  
Contact: sreekar2858.tech  
GitHub: [github.com/sreekar2858](https://github.com/sreekar2858)  
Portfolio: [sreekar2858.tech](https://sreekar2858.tech/)

## Features

- **Configurable**: All simulation parameters are set in a single YAML config file.
- **Modular**: Mesh generation, boundary marking, and solver logic are cleanly separated.
- **Exportable**: Meshes, boundary conditions, and mode shapes can be exported for visualization.
- **Extensible**: Easily adapt to new geometries, materials, or solver settings.

## Project Structure

```
dolfinx/
│
├── config/
│   ├── config_boxMesh.yaml  # Configuration for box mesh simulations
│   └── config_gmsh.yaml     # Configuration for GMSH-based simulations
│
├── data/                   # Geometry files (STL, STEP, etc.)
│   ├── cantilever.stl      # Simple cantilever beam
│   ├── cantilever_long.stl # Long cantilever beam
│   ├── cantilever_comet.stl # Comet-shaped cantilever
│   └── cantilever_comet_scaled_0.001.stl # Scaled comet cantilever
│
├── src/
│   ├── __init__.py        # Package initialization
│   └── solver.py          # Core solver logic
│
├── utils/
│   ├── __init__.py        # Package initialization 
│   ├── mesh_utils.py      # Mesh generation and boundary utilities
│   ├── file_utils.py      # Config and file loading
│   └── bc.py              # Boundary condition helpers
│
├── tests/
│   ├── test.ipynb         # Example/test notebook
│   └── cantilever.py      # Test script for cantilever analysis
│
├── output/                # Generated output files
│   ├── *.vtu              # VTK unstructured grid files
│   ├── *.pvtu             # Parallel VTK unstructured grid files
│   ├── *.h5               # HDF5 data files
│   └── *.xdmf             # XDMF files for visualization
│
├── __init__.py            # Package initialization
├── main.py                # Main entry point
├── setup.py               # Package setup script
├── requirements.txt       # Package dependencies
├── LICENSE                # GPL-3.0 license
├── CONTRIBUTING.md        # Contributing guidelines
└── README.md              # This documentation file
```

## Quick Start

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/sreekar2858/dolfinx.git
   cd dolfinx
   ```

2. **Install DOLFINx and dependencies**
   
   DOLFINx is not available on PyPI and must be installed using conda:
   
   ```bash
   # Create a new conda environment
   conda create -n fenicsx-env
   
   # Activate the environment
   conda activate fenicsx-env
   
   # Install DOLFINx and dependencies
   conda install -c conda-forge fenics-dolfinx mpich pyvista
   
   # Install remaining Python dependencies
   pip install -r requirements.txt
   ```
   
   For more detailed installation instructions or alternative installation methods,
   please visit the [FEniCS Project download page](https://fenicsproject.org/download/).

### Running Simulations

1. **Choose a configuration file**
   - `config_boxMesh.yaml`: For simple box mesh simulations
   - `config_gmsh.yaml`: For STL/STEP geometry-based simulations

2. **Run the simulation**
   ```bash
   # Using the default configuration
   python main.py
   
   # Or specify a config file
   python main.py -c config/config_boxMesh.yaml
   ```

3. **View the results**
   - Output files are saved to the `output/` directory
   - Mesh files (`.vtu`): View in ParaView or other VTK viewers
   - Mode shapes (`.xdmf`): View in ParaView with the accompanying `.h5` files

### Example Output

For each simulation, the following files are typically generated:
- Mesh representation (e.g., `cantilever_box_mesh.vtu`)
- Boundary condition markers (e.g., `cantilever_box_bc_marked.vtu`)
- Mode shapes for each computed eigenmode (e.g., `cantilever_box_modes_1.xdmf`, etc.)

## Configuration

The solver supports two types of configurations in the `config/` directory:

### Box Mesh Configuration (`config_boxMesh.yaml`)

For simple parametric box meshes:

```yaml
material:
  E: 1e5          # Young's modulus [Pa]
  nu: 0.0         # Poisson's ratio
  density: 1e-3   # Density [kg/m³]

geometry:
  dimensions:
    x: [-0.25, 0.25]  # X-dimension bounds [min, max]
    y: [0, 20]        # Y-dimension bounds [min, max]
    z: [-0.5, 0.5]    # Z-dimension bounds [min, max]
  use_boxmesh: true   # Use built-in box mesh generator

mesh:
  divisions:
    Nx: 5            # Number of elements in X direction
    Ny: 100          # Number of elements in Y direction
    Nz: 10           # Number of elements in Z direction
  element_order: 1    # Order of finite elements

solver:
  name: "cantilever"  # Name prefix for output files
  num_modes: 6        # Number of eigenvalues/modes to compute
  tol: 1e-8           # Tolerance for eigenvalue solver
  max_it: 100         # Maximum iterations for eigenvalue solver

exports:
  mesh:
    filename: "cantilever_box_mesh.vtu"       # Output mesh file
  bc_marker:
    filename: "cantilever_box_bc_marked.vtu"  # Output boundary condition file
  modes:
    filename: "cantilever_box_modes"          # Prefix for mode shape files
```

### GMSH Configuration (`config_gmsh.yaml`)

For complex geometries from STL/STEP files:

```yaml
material:
  E: 1e5          # Young's modulus [Pa]
  nu: 0.0         # Poisson's ratio
  density: 1e-3   # Density [kg/m³]

geometry:
  input_file: "data/cantilever_comet.stl"  # Path to geometry file
  scale: 1e-3                             # Scale factor for geometry
  use_boxmesh: false                      # Use GMSH instead of box mesh

mesh:
  min_size: 0.1                           # Minimum mesh element size
  max_size: 0.1                           # Maximum mesh element size
  element_order: 1                        # Order of finite elements

solver:
  name: "cantilever"                      # Name prefix for output files
  num_modes: 6                            # Number of eigenvalues/modes to compute
  tol: 1e-8                               # Tolerance for eigenvalue solver
  max_it: 100                             # Maximum iterations for eigenvalue solver

exports:
  mesh:
    filename: "cantilever_gmsh_mesh.vtu"       # Output mesh file
  bc_marker:
    filename: "cantilever_gmsh_bc_marked.vtu"  # Output boundary condition file
  modes:
    filename: "cantilever_gmsh_modes"          # Prefix for mode shape files
```

## Codebase Structure and Execution

### Key Components

1. **Configuration System**
   - All simulation parameters are defined in YAML files under `config/`
   - Loaded and parsed by `utils/file_utils.py`

2. **Mesh Handling**
   - `utils/mesh_utils.py`: Creates meshes from geometries or box meshes
   - Supports both analytical box meshes and STL/STEP imports via GMSH
   - Handles boundary marking and facet identification

3. **Physics Setup**
   - `src/solver.py`: Defines the elasticity problem
   - Sets up material properties, boundary conditions
   - Creates function spaces and forms

4. **Modal Analysis**
   - Eigenvalue solver configuration in `src/solver.py`
   - Uses SLEPc for efficient eigenvalue extraction
   - Computes natural frequencies and mode shapes

5. **Output and Visualization**
   - Results exported as VTK, XDMF, and HDF5 files
   - Compatible with ParaView and other scientific visualization tools

### Execution Flow

1. `main.py` parses command-line arguments and loads the configuration
2. `run_modal_analysis()` from `src/solver.py` is called with the config path
3. The mesh is created based on configuration (box mesh or STL)
4. Boundary conditions are applied to the fixed end of the cantilever
5. The stiffness matrix (A) and mass matrix (M) are assembled
6. The eigenvalue problem Ax = λMx is solved using SLEPc
7. Natural frequencies and mode shapes are extracted and exported

## Extending

- To use a different geometry, update the `input_file` in the config.
- To change boundary conditions, modify the relevant function in `utils/mesh_utils.py` and reference it in `src/solver.py`.
- To add new exports, update the `exports` section in the config and the corresponding logic in `src/solver.py`.
- To implement different material models, modify the stress-strain relationship in `src/solver.py`.

## Testing and Development

### Test Scripts

- `tests/cantilever.py`: A standalone script that performs modal analysis on a cantilever beam
  - Can be used to verify the solver against analytical solutions
  - Demonstrates how to use the DOLFINx API directly

To run the test script:

```bash
python tests/cantilever.py
```

### Visualization

The generated mode shapes can be visualized using ParaView:

1. Open ParaView
2. Load the `.vtu` files from the `output/` directory
3. Apply the "Warp By Vector" filter to visualize the mode shapes with exaggerated displacement
