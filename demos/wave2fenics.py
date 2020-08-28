import time
import dolfin
from raschii import get_wave_model, check_breaking_criteria


def wave2fenics(wave, t, rangex, rangey, Nx, Ny, xdmf=False, plot=False):
    """
    Show how to construct a FEniCS velocity field based on the wave velocity
    field. Tested with FEniCS 2018.1 (pre-release dev version)
    """
    t1 = time.time()

    start = dolfin.Point(rangex[0], rangey[0])
    end = dolfin.Point(rangex[1], rangey[1])
    mesh = dolfin.RectangleMesh(dolfin.MPI.comm_world, start, end, Nx, Ny)

    # Velocity and divergence
    u0, u1, div, tot_flux = make_vel_div(wave, mesh, t)

    print("DONE in %.2f seconds" % (time.time() - t1))

    if xdmf:
        xdmf_file = "divergence.xdmf"
        print("    Writing XDFM file %r" % xdmf_file)

        # Add a colour function to the XDMF file to show the wave profile
        c = make_colour_function(wave, mesh, t)

        # Write the results to XDMF
        with dolfin.XDMFFile(mesh.mpi_comm(), xdmf_file) as xdmf:
            xdmf.parameters["rewrite_function_mesh"] = False
            xdmf.parameters["functions_share_mesh"] = True
            u0.rename("u0", "u0")
            xdmf.write(u0, 0.0)
            u1.rename("u1", "u1")
            xdmf.write(u1, 0.0)
            div.rename("div", "div")
            xdmf.write(div, 0.0)
            c.rename("c", "c")
            xdmf.write(c, 0.0)

    print("\nStatistics of the results:")
    print("    Total flux across the domain boundaries in FEM is    % .4e" % tot_flux)
    print(
        "    Maximum absolute divergence of the velocity in FEM is %.4e"
        % abs(div.vector().get_local()).max()
    )

    # Compare FEniCS interpolated values with the Raschii Python values
    V = u0.function_space()
    dofs_x = V.tabulate_dof_coordinates().reshape((-1, 2))
    vel_py = wave.velocity(dofs_x[:, 0], dofs_x[:, 1], t)
    vel_f0 = u0.vector().get_local()
    vel_f1 = u1.vector().get_local()
    print(
        "    Maximum absolute difference at FEniCS dofs of u0 is   %.4e"
        % abs(vel_py[:, 0] - vel_f0).max()
    )
    print(
        "    Maximum absolute difference at FEniCS dofs of u1 is   %.4e"
        % abs(vel_py[:, 1] - vel_f1).max()
    )

    if not plot:
        return

    png_file = "divergence.png"
    print("\nPlotting with matplotlib to", png_file)

    import numpy
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot

    # Plot the FEniCS interpolated velocity divergence
    eps = 1e-6
    xpos = numpy.linspace(rangex[0] + eps * 2, rangex[1] - eps * 2, Nx * 3)
    ypos = numpy.linspace(rangey[0] + eps * 2, rangey[1] - eps * 2, Ny * 3)
    X, Y = numpy.meshgrid(xpos, ypos)
    xr, yr = X.ravel(), Y.ravel()
    divF = numpy.zeros_like(xr)
    pos = numpy.zeros(2, float)
    val = numpy.zeros(1, float)
    for i in range(xr.size):
        pos[:] = (xr[i], yr[i])
        div.eval(val, pos)
        divF[i] = val[0]

    eta = wave.surface_elevation(xpos, t)

    fig = pyplot.figure(figsize=(10, 6))

    ax = fig.add_subplot(111)
    c = ax.contourf(X, Y, abs(divF).reshape(X.shape))
    fig.colorbar(c, ax=ax)
    ax.plot(xpos, eta)

    fig.tight_layout()
    fig.savefig(png_file)


def make_vel_div(wave, mesh, t):
    """
    Interpolate velocities into DG a Lagrange function space and project the
    divergence of the resulting velocity vector into the same space
    """
    V = dolfin.FunctionSpace(mesh, "DG", 2)
    u0 = dolfin.Function(V)
    u1 = dolfin.Function(V)
    vel = dolfin.as_vector([u0, u1])

    # Get C++ expressions and convert to x-y plane from x-z plane
    cpp_x, cpp_y = wave.velocity_cpp()
    cpp_x = cpp_x.replace("x[2]", "x[1]")
    cpp_y = cpp_y.replace("x[2]", "x[1]")

    # Interpolate expressions
    e0 = dolfin.Expression(cpp_x, t=t, degree=2)
    e1 = dolfin.Expression(cpp_y, t=t, degree=2)
    u0.interpolate(e0)
    u1.interpolate(e1)

    # Project divergence
    u, v = dolfin.TrialFunction(V), dolfin.TestFunction(V)
    div = dolfin.Function(V)
    A = dolfin.assemble(u * v * dolfin.dx)
    b = dolfin.assemble(dolfin.div(vel) * v * dolfin.dx)
    dolfin.solve(A, div.vector(), b)

    # Total flux across the boundaries
    n = dolfin.FacetNormal(mesh)
    tot_flux = dolfin.assemble(dolfin.dot(vel, n) * dolfin.ds)

    return u0, u1, div, tot_flux


def make_colour_function(wave, mesh, t):
    """
    Make colour function and project to DG0 (piecewise constants)
    """
    V0 = dolfin.FunctionSpace(mesh, "DG", 0)
    quad_degree = 6
    element = dolfin.FiniteElement(
        "Quadrature", mesh.ufl_cell(), quad_degree, quad_scheme="default"
    )
    ec = dolfin.Expression("x[1] < (%s) ? 1.0 : 0.0" % wave.elevation_cpp(), element=element, t=t)
    u0, v0 = dolfin.TrialFunction(V0), dolfin.TestFunction(V0)
    c = dolfin.Function(V0)
    A0 = dolfin.assemble(u0 * v0 * dolfin.dx)
    b0 = dolfin.assemble(ec * v0 * dolfin.dx(metadata={"quadrature_degree": quad_degree}))
    dolfin.solve(A0, c.vector(), b0)
    return c


def main():
    """
    Parse command line arguments and produce XDMF file of the velocity field,
    the velocity field divergence and the colour function using FEniCS DOLFIN
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="raschii_wave2fenics", description="Plot the divergence of the velocity"
    )
    parser.add_argument("wave_type", help='Name of the wave model, e.g., "Fenton"')
    parser.add_argument("wave_height", help="Wave height", type=float)
    parser.add_argument("water_depth", help="The still water depth", type=float)
    parser.add_argument("wave_length", help="Distance between peaks", type=float)
    parser.add_argument("-N", type=int, default=10, help="Approximation order")
    parser.add_argument("--Nx", type=int, default=101, help="#cells in x-direction")
    parser.add_argument("--Ny", type=int, default=101, help="#cells in y-direction")
    parser.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="Allow exceeding breaking criteria",
    )
    parser.add_argument(
        "-a",
        "--depth_air",
        default=0.0,
        type=float,
        help="Include velocities in the air phase if this is "
        'greater than 0 (distance to "lid" above the air).',
    )
    parser.add_argument(
        "-b",
        "--blend_distance",
        default=-1,
        type=float,
        help="Blend the water and air stream functions a "
        "distance up to improve continuity of velocities",
    )
    parser.add_argument("-t", "--time", default=0.0, type=float, help="The time instance to plot")
    parser.add_argument("--ymin", default=None, type=float, help="Lower vertical axis limit")
    parser.add_argument("--ymax", default=None, type=float, help="Upper vertical axis limit")
    parser.add_argument(
        "--plot",
        default=False,
        action="store_true",
        help="Show matplotlib plot of eval-ed dolfin.Functions",
    )
    parser.add_argument("--xdmf", default=False, action="store_true", help="Produce xdmf file")
    args = parser.parse_args()

    err, warn = check_breaking_criteria(args.wave_height, args.water_depth, args.wave_length)
    if err:
        print(err)
    if warn:
        print(warn)
    if err and not args.force:
        exit(1)

    t1 = time.time()
    print("\nGenerating wave, may take some time ...")
    WaveClass, AirClass = get_wave_model(args.wave_type)
    wave_args = dict(height=args.wave_height, depth=args.water_depth, length=args.wave_length)

    if "N" in WaveClass.required_input:
        wave_args["N"] = args.N
    if "air" in WaveClass.optional_input and AirClass is not None:
        blend_distance = None if args.blend_distance < 0 else args.blend_distance
        wave_args["air"] = AirClass(args.depth_air, blend_distance)

    # Print wave input, then generate the wave
    for a in sorted(wave_args):
        print("%13s: %5r" % (a, wave_args[a]))
    wave = WaveClass(**wave_args)
    print("DONE in %.2f seconds" % (time.time() - t1))

    # Determine the dolfin mesh size
    height = args.water_depth + max(args.depth_air, args.wave_height * 2)
    length = args.wave_length * 2
    xmin = 0
    xmax = length
    ymin = 0 if args.ymin is None else args.ymin
    ymax = height if args.ymax is None else args.ymax
    if args.ymax is None:
        ymax = args.water_depth + max(args.depth_air, args.wave_height * 2)

    print("\nConverting to FEniCS C++ expressions and interpolating to FEM ...")
    wave2fenics(
        wave,
        args.time,
        (xmin, xmax),
        (ymin, ymax),
        args.Nx,
        args.Ny,
        plot=args.plot,
        xdmf=args.xdmf,
    )


if __name__ == "__main__":
    print("RUNNING wave2fenics:")
    main()
    print("\nDONE")
