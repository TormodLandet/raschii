import numpy
from matplotlib import pyplot
from raschii import get_wave_model, WAVE_MODELS, check_breaking_criteria


def plot_wave(wave_type, height, depth, length, N, linear):
    """
    Plot the wave with the given parameters
    """
    WaveClass = get_wave_model(wave_type)
    args = dict(height=height, depth=depth, length=length)
    if 'N' in WaveClass.required_input:
        args['N'] = N
    wave = WaveClass(**args)
    
    x = numpy.linspace(0, length/2, 200)
    eta = wave.surface_elevation(x)
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    
    # Plot the linear surface elevation
    if linear:
        eta_lin = depth + height/2 * numpy.cos(2 * numpy.pi * x / length)
        ax.plot(x, eta_lin, 'g:')
    
    # Plot surface elevation (with colocation points)
    ax.plot(x, eta)
    ax.plot(wave.x, wave.eta, 'kx', ms=2)
    
    ymax = depth + height
    ymin = max(depth - length * 2, 0)    # Include only 2 wave lengths dept,
    ymin = max(depth - height *10, ymin) # but no more than 10 wave heights
    ax.set_ylim(ymin, ymax)
    
    pyplot.show()


if __name__ == '__main__':
    # Get command line arguments
    import argparse
    parser = argparse.ArgumentParser(prog='raschii_plot_wave',
                                     description='Plot a wave')
    parser.add_argument('wave_type', help='Name of the wave model',
                        choices=WAVE_MODELS.keys())
    parser.add_argument('wave_height', help='Wave height', type=float)
    parser.add_argument('water_depth', help='The still water depth', type=float)
    parser.add_argument('wave_length', help='Distance between peaks', type=float)
    parser.add_argument('-N', type=int, default=10, help='Approximation order')
    parser.add_argument('-l', '--linear', default=False, action='store_true',
                        help='Plot the linear wave elevation')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Allow exceeding breaking criteria')
    args = parser.parse_args()
    
    err, warn = check_breaking_criteria(args.wave_height, args.water_depth,
                                        args.wave_length)
    if err: print(err)
    if warn: print(warn)
    if err and not args.force:
        exit(1)
    
    plot_wave(args.wave_type, args.wave_height, args.water_depth,
              args.wave_length, args.N, args.linear)
