"""Run the single-pad baseline and numerical consistency checks (not ANSYS parity)."""
import argparse
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from dolfinx.mesh import create_rectangle, CellType
from thermoelastic_solver import ReferenceGeometry, cu_at_250C, solve


def create_baseline_mesh(g, nr=40, nz=40):
    if nr < 4 or nr % 2 or nz < 4:
        raise ValueError('nr must be even and >=4; nz must be >=4')
    # Coordinate mapping preserves mesh rows at all material interfaces.
    gap_cells = 2 if g.recess_nm > 0 else 0
    total_z = 2*nz + gap_cells
    mesh = create_rectangle(MPI.COMM_WORLD, [[0., 0.], [1., 1.]],
                            [nr, total_z], CellType.triangle)
    x = mesh.geometry.x
    x[:,0] = np.interp(x[:,0], [0., .5, 1.],
                       [0., g.pad_radius_nm/g.L0_nm, g.outer_radius_nm/g.L0_nm])
    knots = [0., nz/total_z, 2*nz/total_z]
    heights = [0., g.base_thickness_nm/g.L0_nm, g.cu_top]
    if gap_cells:
        knots.append(1.)
        heights.append(g.top)
    x[:,1] = np.interp(x[:,1], knots, heights)
    return mesh


def run(g, mesh, temperature, residual=1e3):
    return solve(None, mesh, cu_at_250C(), temperature, g.L0_m,
                 g.pad_radius_nm/g.L0_nm, (g.pad_height_nm-g.recess_nm)/g.L0_nm,
                 geometry=g, residual_void_E_Pa=residual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nr', type=int, default=40)
    parser.add_argument('--nz', type=int, default=40)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--ansys', type=Path, help='Combined ANSYS CSV in physical reference coordinates')
    parser.add_argument('--output', type=Path, default=Path('baseline-results'))
    args = parser.parse_args()
    g = ReferenceGeometry()
    mesh = create_baseline_mesh(g, args.nr, args.nz)
    result = run(g, mesh, 523.15)
    if args.check:
        zero = run(g, mesh, 293.15)
        half = run(g, mesh, 408.15)
        assert np.allclose(zero.u_dimless_field.x.array, 0., atol=1e-14)
        assert np.allclose(half.u_dimless_field.x.array, result.u_dimless_field.x.array/2,
                           rtol=1e-8, atol=1e-12)
        assert np.allclose(half.W_field.x.array, result.W_field.x.array/4, rtol=1e-7, atol=1e-5)
        coords = result.u_dimless_field.function_space.tabulate_dof_coordinates()
        values = result.u_dimless_field.x.array.reshape(-1,2)
        radial_fixed = np.isclose(coords[:,0], 0.) | np.isclose(coords[:,0], 1.) | np.isclose(coords[:,1], 0.)
        assert np.allclose(values[radial_fixed,0], 0., atol=1e-12)
        assert np.allclose(values[np.isclose(coords[:,1],0.),1], 0., atol=1e-12)
        assert np.all(np.isfinite(result.W_full_field.x.array))
        assert np.min(result.W_full_field.x.array) >= -1e-6
        # The physical problem must be independent of the normalization length.
        from dataclasses import replace
        scaled_g = replace(g, L0_nm=1000.)
        scaled = run(scaled_g, create_baseline_mesh(scaled_g, args.nr, args.nz), 523.15)
        assert np.allclose(scaled.cu_surface_m, result.cu_surface_m, rtol=1e-8, atol=1e-15)
        # Homogeneous free expansion has the exact solution u = alpha*dT*X.
        from dolfinx import fem
        from thermoelastic_solver import MechanicalBCs
        phase = fem.Function(fem.functionspace(mesh, ("Lagrange", 1)))
        phase.x.array[:] = 1.
        material = cu_at_250C()
        free = solve(phase, mesh, material, 523.15, g.L0_m,
                     g.pad_radius_nm/g.L0_nm, (g.pad_height_nm-g.recess_nm)/g.L0_nm,
                     geometry=g, material_dielectric=material,
                     boundary_conditions=MechanicalBCs(bottom_ur_fixed=False, outer_ur_fixed=False))
        expected = material.alpha_1_K*230.*coords[:,:2]
        assert np.allclose(free.u_dimless_field.x.array.reshape(-1,2), expected, atol=1e-10)
        assert np.max(np.abs(free.W_full_field.x.array)) < 1e-6
        if mesh.comm.rank == 0:
            print('PASS: zero load, thermal scaling, energy scaling, BCs, finite energy, length normalization, analytic free expansion')
    if mesh.comm.rank == 0:
        args.output.mkdir(parents=True, exist_ok=True)
        if args.ansys:
            import csv
            with args.ansys.open(newline='') as stream:
                records = list(csv.DictReader(stream))
            reference = np.array([[float(row['X Location (m)']), float(row['Y Location (m)']),
                                   float(row['Directional Deformation Y (m)'])] for row in records])
            mask = ((reference[:,0] <= g.pad_radius_nm*1e-9+1e-12)
                    & (np.abs(reference[:,1]-g.cu_top*g.L0_m) < 1e-12))
            reference = reference[mask]
            if not len(reference):
                raise ValueError('No ANSYS nodes on the reference Cu top surface')
            radii, inverse = np.unique(reference[:,0], return_inverse=True)
            displacement = np.bincount(inverse, weights=reference[:,2])/np.bincount(inverse)
            prediction = np.interp(radii, result.UY_r_dimless*g.L0_m, result.UY_dimless*g.L0_m)
            error = prediction-displacement
            result.stats['ansys_surface_comparison'] = {
                'source': str(args.ansys), 'n_points': len(radii),
                'center_uz_nm': float(displacement[0]*1e9),
                'edge_uz_nm': float(displacement[-1]*1e9),
                'rmse_nm': float(np.sqrt(np.mean(error**2))*1e9),
                'max_abs_error_nm': float(np.max(np.abs(error))*1e9)}
            np.savetxt(args.output/'ansys_surface_comparison.csv',
                       np.column_stack((radii, displacement, prediction, error)), delimiter=',',
                       header='r_m,ansys_uz_m,fenicsx_uz_m,error_m', comments='')
        (args.output/'summary.json').write_text(json.dumps(result.stats, indent=2)+'\n')
        np.savetxt(args.output/'cu_surface.csv', np.column_stack((result.UY_r_dimless*g.L0_m,
                   result.UY_dimless*g.L0_m, result.cu_surface_m)), delimiter=',',
                   header='r_reference_m,uz_m,r_deformed_m,z_deformed_m', comments='')
        np.savetxt(args.output/'dielectric_surface.csv', result.dielectric_surface_m,
                   delimiter=',', header='r_deformed_m,z_deformed_m', comments='')
        print(json.dumps(result.stats, indent=2))
        print('ANSYS parity remains unverified.')


if __name__ == '__main__':
    main()
