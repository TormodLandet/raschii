import numpy
from math import pi, sinh, cosh, tanh
from matplotlib import pyplot
from raschii import get_wave_model, check_breaking_criteria


def plot_wave(model_names, height, depth, length, N):
    """
    Plot waves with the given parameters
    """
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    # Print some summary info about the waves
    head_format = '%15s  %10s %10s %10s %10s'
    info_format = '%15s  %10.3e %10.3e %10.3e %10.3e'
    print(head_format % ('', 'c', 'eta_max', 'eta_min', 'vel_max'))
    
    for model_name in model_names:
        WaveClass = get_wave_model(model_name)
        args = dict(height=height, depth=depth, length=length)
        if 'N' in WaveClass.required_input:
            args['N'] = N
        wave = WaveClass(**args)
        
        # Vertical axis limits
        ymax = depth + height * 1.5
        ymin = max(depth - length * 2, 0)    # Include only 2 wave lengths dept,
        ymin = max(depth - height * 10, ymin) # but no more than 10 wave heights
        
        # Get elevation
        x = numpy.linspace(0, length/2, 200)
        eta = wave.surface_elevation(x)
        
        # Get velocity
        X, Y = numpy.meshgrid(numpy.linspace(x[0], x[-1], 11),
                              numpy.linspace(ymin, eta.max(), 11))
        vel = wave.velocity(X.ravel(), Y.ravel())
        U = vel[:,0].reshape(X.shape)
        V = vel[:,1].reshape(X.shape)
        
        # Crest velocity scale (Airy value)
        k_scale = 2 * pi / length
        g = 9.81
        w_scale = (k_scale * g * tanh(k_scale * depth))**0.5
        Uscale = w_scale * height / 2 * (cosh(k_scale * (depth + height / 2)) /
                                         sinh(k_scale * depth))
        
        # Plot surface elevation
        eta_line, = ax.plot(x, eta, label=model_name)
        eta_color = eta_line.get_color()
        
        # Plot colocation points
        if hasattr(wave, 'eta'):
            ax.plot(wave.x, wave.eta, 'kx', ms=2)
        
        # Plot the velocity components
        ax.quiver(X, Y, U, V, scale_units='height', scale=Uscale*8,
                  color=eta_color, width=3.5e-3, label=None)
        
        # Print some info
        vel_max = ((vel[:,0]**2 + vel[:,1]**2)**0.5).max()
        print(info_format % (model_name, wave.c, eta.max(), eta.min(), vel_max))
    
    # Print the velocity scale
    print(head_format % ('scale', '', '', '', '%10.3e' % Uscale))
    
    # Show a legend if there are more than one wave model being plotted
    if len(model_names) > 1:    
        ax.legend(loc='lower right')
        ax.set_title('Wave models %s' % ' and '.join(model_names))
    else:
        ax.set_title('Wave model %s' % model_name)
    
    ax.set_ylim(ymin, ymax)
    pyplot.show()


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
    args = parser.parse_args()
    
    err, warn = check_breaking_criteria(args.wave_height, args.water_depth,
                                        args.wave_length)
    if err: print(err)
    if warn: print(warn)
    if err and not args.force:
        exit(1)
    
    model_names = args.wave_type.split(',')
    plot_wave(model_names, args.wave_height, args.water_depth,
              args.wave_length, args.N)
