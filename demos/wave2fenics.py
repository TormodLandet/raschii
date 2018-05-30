import dolfin
from raschii import get_wave_model, check_breaking_criteria


def wave2fenics(wave, rangex, rangey, Nx=101, Ny=51):
    """
    Show how to construct a FEniCS velocity field based on the wave velocity
    field. Tested with FEniCS 2018.1 (pre-release dev version)
    """
    start = dolfin.Point(rangex[0], rangey[0])
    end = dolfin.Point(rangex[1], rangey[1])
    mesh = dolfin.RectangleMesh(dolfin.MPI.comm_world, start, end, Nx, Ny)
    
    # Velocity and divergence
    u0, u1, div = make_vel_div(wave, mesh)
    
    # Colour function
    c = make_colour_function(wave, mesh)
    
    # Plot results
    file_name = 'divergence.xdmf'
    print('    Writing XDFM file %r' % file_name)
    with dolfin.XDMFFile(mesh.mpi_comm(), file_name) as xdmf:
        xdmf.parameters['rewrite_function_mesh'] = False
        xdmf.parameters['functions_share_mesh'] = True
        u0.rename('u0', 'u0')
        xdmf.write(u0, 0.0)
        u1.rename('u1', 'u1')
        xdmf.write(u1, 0.0)
        div.rename('div', 'div')
        xdmf.write(div, 0.0)
        c.rename('c', 'c')
        xdmf.write(c, 0.0)


def make_vel_div(wave, mesh):
    """
    Interpolate velocities into DG a Lagrange function space and project the
    divergence of the resulting velocity vector into the same space
    """
    V = dolfin.FunctionSpace(mesh, 'DG', 2)
    u0 = dolfin.Function(V)
    u1 = dolfin.Function(V)
    vel = dolfin.as_vector([u0, u1])
    
    # Get C++ expressions and convert to x-y plane from x-z plane
    cpp_x, cpp_y = wave.velocity_cpp()
    cpp_x = cpp_x.replace('x[2]', 'x[1]')
    cpp_y = cpp_y.replace('x[2]', 'x[1]')
    
    # Interpolate expressions
    e0 = dolfin.Expression(cpp_x, t=0, degree=2)
    e1 = dolfin.Expression(cpp_y, t=0, degree=2)
    u0.interpolate(e0)
    u1.interpolate(e1)
    
    # Project divergence
    u, v = dolfin.TrialFunction(V), dolfin.TestFunction(V)
    div = dolfin.Function(V)
    A = dolfin.assemble(u * v * dolfin.dx)
    b = dolfin.assemble(dolfin.div(vel) * v * dolfin.dx)
    dolfin.solve(A, div.vector(), b)
    
    return u0, u1, div


def make_colour_function(wave, mesh):
    """
    Make colour function and project to DG0 (piecewise constants)
    """
    V0 = dolfin.FunctionSpace(mesh, 'DG', 0)
    quad_degree = 6
    element = dolfin.FiniteElement('Quadrature', mesh.ufl_cell(),
                                       quad_degree, quad_scheme='default')
    ec = dolfin.Expression('x[1] < (%s) ? 1.0 : 0.0' % wave.elevation_cpp(),
                           element=element, t=0)
    u0, v0 = dolfin.TrialFunction(V0), dolfin.TestFunction(V0)
    c = dolfin.Function(V0)
    A0 = dolfin.assemble(u0 * v0 * dolfin.dx)
    b0 = dolfin.assemble(ec * v0 * dolfin.dx(metadata={
        'quadrature_degree': quad_degree}))
    dolfin.solve(A0, c.vector(), b0)
    return c


def main():
    """
    Parse command line arguments and produce XDMF file of the velocity field,
    the velocity field divergence and the colour function using FEniCS DOLFIN
    """
    import argparse
    parser = argparse.ArgumentParser(prog='raschii_wave2fenics',
                                     description='Plot the divergence of the velocity')
    parser.add_argument('wave_type', help='Name of the wave model, e.g., "Fenton"')
    parser.add_argument('wave_height', help='Wave height', type=float)
    parser.add_argument('water_depth', help='The still water depth', type=float)
    parser.add_argument('wave_length', help='Distance between peaks', type=float)
    parser.add_argument('-N', type=int, default=10, help='Approximation order')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Allow exceeding breaking criteria')
    parser.add_argument('-a', '--depth_air', default=0.0, type=float,
                        help='Include velocities in the air phase if this is '
                             'greater than 0 (distance to "lid" above the air).')
    parser.add_argument('-b', '--blend_distance', default=-1, type=float,
                        help='Blend the water and air stream functions a '
                             'distance up to improve continuity of velocities')
    parser.add_argument('-t', '--time', default=0.0, type=float,
                        help='The time instance to plot')
    parser.add_argument('--ymin', default=None, type=float,
                        help='Lower vertical axis limit')
    parser.add_argument('--ymax', default=None, type=float,
                        help='Upper vertical axis limit')
    args = parser.parse_args()
    
    err, warn = check_breaking_criteria(args.wave_height, args.water_depth,
                                        args.wave_length)
    if err: print(err)
    if warn: print(warn)
    if err and not args.force:
        exit(1)
    
    print('\nGenerating wave, may take some time ...')
    WaveClass, AirClass = get_wave_model(args.wave_type)
    wave_args = dict(height=args.wave_height,
                     depth=args.water_depth,
                     length=args.wave_length)
    if 'N' in WaveClass.required_input:
        wave_args['N'] = args.N
    if 'air' in WaveClass.optional_input and AirClass is not None:
        blend_distance = None if args.blend_distance < 0 else args.blend_distance
        wave_args['air'] = AirClass(args.depth_air, blend_distance)
    
    for a in sorted(wave_args):
        print('%13s: %5r' % (a, wave_args[a]))
    
    wave = WaveClass(**wave_args)
    
    print('\nConverting to FEniCS and writing XDMF ...')
    height = args.water_depth + max(args.depth_air, args.wave_height * 2)
    length = args.wave_length * 2 
    xmin = 0
    xmax = length
    ymin = 0 if args.ymin is None else args.ymin
    ymax = height if args.ymax is None else args.ymax    
    if args.ymax is None:
        ymax = args.water_depth + max(args.depth_air, args.wave_height * 2)
    wave2fenics(wave, (xmin, xmax), (ymin, ymax))
    

if __name__ == '__main__':
    print('RUNNING wave2fenics:')
    main()
    print('\nDONE')
