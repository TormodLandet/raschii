import numpy
from math import pi, sinh, cosh, tanh
from matplotlib import pyplot
from raschii import get_wave_model, check_breaking_criteria


def plot_wave(model_names, height, depth, length, N, depth_air, t, Nx=21, Ny=21):
    """
    Plot waves with the given parameters
    """
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    # Print some summary info about the waves
    head_format = '%15s  %10s %10s %10s %10s %10s'
    info_format = '%15s  %10.3e %10.3e %10.3e %10.3e %10.3e'
    print(head_format % ('', 'c', 'T', 'eta_max', 'eta_min', 'vel_max'))
    
    ymax2 = 0
    for model_name in model_names:
        WaveClass = get_wave_model(model_name)
        args = dict(height=height, depth=depth, length=length)
        if 'N' in WaveClass.required_input:
            args['N'] = N
        if 'depth_air' in WaveClass.optional_input and depth_air > 0:
            args['depth_air'] = depth_air
        wave = WaveClass(**args)
        
        # Vertical axis limits
        ymax = depth + height * 1.5
        ymin = max(depth - length * 2, 0)    # Include only 2 wave lengths dept,
        ymin = max(depth - height * 10, ymin) # but no more than 10 wave heights
        
        # Get elevation
        x = numpy.linspace(-length/2, length/2, 200)
        eta = wave.surface_elevation(x, t)
        
        # Get velocity
        xvec = numpy.linspace(-x[-1], x[-1], Nx)
        evec = wave.surface_elevation(xvec, t)
        X, Y = makeXY(xvec, xvec * 0 + ymin, evec, Ny)
        vel = wave.velocity(X.ravel(), Y.ravel(), t)
        U = vel[:,0].reshape(X.shape)
        V = vel[:,1].reshape(X.shape)
        
        # Crest velocity scale (Airy value)
        k_scale = 2 * pi / length
        g = 9.81
        w_scale = (k_scale * g * tanh(k_scale * depth))**0.5
        Uscale = w_scale * height / 2 * (cosh(k_scale * (depth + height / 2)) /
                                         sinh(k_scale * depth))
        scale = Uscale * 8
        
        # Plot surface elevation
        eta_line, = ax.plot(x, eta, label=model_name, zorder=1000, lw=2)
        eta_color = eta_line.get_color()
        
        # Plot colocation points
        if hasattr(wave, 'eta') and t == 0:
            ax.plot(wave.x, wave.eta, 'kx', ms=2)
        
        # Plot the velocity components
        ax.quiver(X, Y, U, V, scale_units='height', scale=scale,
                  color=eta_color, width=3.5e-3, label=None, alpha=0.6)
        
        # Plot the velocities in the air
        if hasattr(wave, 'air') and wave.air is not None:
            ymax2 = min(depth + length * 2, depth + depth_air)
            ymax2 = min(depth + height * 10, ymax2)
            X2, Y2 = makeXY(xvec, evec + height / 8, xvec * 0 + ymax2, Ny)
            vel = wave.air.velocity(X2.ravel(), Y2.ravel(), t)
            U2 = vel[:,0].reshape(X.shape)
            V2 = vel[:,1].reshape(X.shape)
            ax.quiver(X2, Y2, U2, V2, scale_units='height', scale=scale,
                      color=eta_color, width=2.0e-3, label=None, alpha=0.6)
        
        # Print some info
        vel_max = ((vel[:,0]**2 + vel[:,1]**2)**0.5).max()
        print(info_format % (model_name, wave.c, 2 * pi / (wave.k * wave.c),
                             eta.max(), eta.min(), vel_max))
    
    # Print the velocity scale
    print(head_format % ('scale', '', '', '', '', '%10.3e' % Uscale))
    
    # Show a legend if there are more than one wave model being plotted
    if len(model_names) > 1:    
        ax.legend(loc='lower right')
        ax.set_title('Wave models %s' % ' and '.join(model_names))
    else:
        ax.set_title('Wave model %s' % model_name)
    
    ax.set_ylim(ymin, max(ymax, ymax2))
    pyplot.show()


def makeXY(x, ymin, ymax, Ny):
    """
    Return X and Y matrices where X positions are spaced according to the input
    x vector and y is linearly interpolated by the boundaries that are given for
    each x position 
    """
    fac = numpy.linspace(0, 1, Ny)
    X, F = numpy.meshgrid(x, fac)
    Y = numpy.zeros_like(X)
    for i in range(x.size):
        for j in range(Ny):
            Y[i,j] = (1 - F[i,j]) * ymin[j] + F[i,j] * ymax[j]
    return X, Y


if __name__ == '__main__':
    # Get command line arguments
    import argparse
    parser = argparse.ArgumentParser(prog='raschii_plot_wave',
                                     description='Plot a wave')
    parser.add_argument('wave_type', help='Name of the wave model')
    parser.add_argument('wave_height', help='Wave height', type=float)
    parser.add_argument('water_depth', help='The still water depth', type=float)
    parser.add_argument('wave_length', help='Distance between peaks', type=float)
    parser.add_argument('-N', type=int, default=10, help='Approximation order')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Allow exceeding breaking criteria')
    parser.add_argument('-a', '--depth_air', default=0.0, type=float,
                        help='Include velocities in the air phase if this is '
                             'greater than 0')
    parser.add_argument('-t', '--time', default=0.0, type=float,
                        help='The time instance to plot')
    args = parser.parse_args()
    
    err, warn = check_breaking_criteria(args.wave_height, args.water_depth,
                                        args.wave_length)
    if err: print(err)
    if warn: print(warn)
    if err and not args.force:
        exit(1)
    
    model_names = args.wave_type.split(',')
    plot_wave(model_names, args.wave_height, args.water_depth,
              args.wave_length, args.N, args.depth_air, args.time)
