import numpy
from math import pi, sinh, cosh, tanh
from matplotlib import pyplot
from raschii import get_wave_model, check_breaking_criteria


def plot_wave(
    model_names,
    height,
    depth,
    length,
    N,
    depth_air,
    blend_height,
    t,
    Nx=21,
    Ny=21,
    plot_quiver=False,
    ymin=None,
    ymax=None,
):
    """
    Plot waves with the given parameters
    """
    # Figure for plot of wave with quiver
    fig1, ax1 = pyplot.subplots(1, 1)

    # Figure for plot of horizontal velocities
    fig2, ax2s = pyplot.subplots(2, 5)
    ax2s = list(ax2s.ravel())

    # Print some summary info about the waves
    head_format = "%20s  %10s %10s %10s %10s %10s %10s %10s"
    info_format = "%20s  %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e"
    print(head_format % ("", "c", "T", "eta_max", "eta_min", "vel_crest", "vel_trough", "vel_max"))

    # Vertical axis limits
    ymin_user, ymax_user = ymin, ymax
    ymax = depth + height * 1.5
    ymin = max(depth - length * 2, 0)  # Include only 2 wave lengths dept,
    ymin = max(depth - height * 10, ymin)  # but no more than 10 wave heights

    # Horizontal velocity axis limits
    vx_max = 0
    vz_max = 0

    # Increase upper vertical limit if plotting the air phase
    if depth_air > 0:
        ymax = min(depth + length * 2, depth + depth_air)
        ymax = min(depth + height * 10, ymax)

    # Let the user override the axis limits
    ymin = ymin_user if ymin_user is not None else ymin
    ymax = ymax_user if ymax_user is not None else ymax

    # Plot each of the specified wave model
    warnings = ""
    for model_name in model_names:
        WaveClass, AirClass = get_wave_model(model_name)
        args = dict(height=height, depth=depth, length=length)
        if "N" in WaveClass.required_input:
            args["N"] = N

        if AirClass is not None and "air" in WaveClass.optional_input:
            blend_height = None if blend_height < 0 else blend_height
            args["air"] = AirClass(depth_air, blend_height)

        wave = WaveClass(**args)
        plot_air = wave.air is not None
        if wave.warnings:
            warnings += "WARNINGS for %s:\n%s" % (model_name, wave.warnings)

        # Get elevation
        assert length > 0
        x = numpy.linspace(-length / 2, length / 2, 200)
        eta = wave.surface_elevation(x, t)

        # Get velocity
        xvec = numpy.linspace(-x[-1], x[-1], Nx)
        evec = wave.surface_elevation(xvec, t)
        X, Y = makeXY(xvec, xvec * 0 + ymin, evec, Ny)
        vel = wave.velocity(X.ravel(), Y.ravel(), t)
        U = vel[:, 0].reshape(X.shape)
        V = vel[:, 1].reshape(X.shape)

        # Crest velocity scale (Airy value)
        k_scale = 2 * pi / length
        g = 9.81
        w_scale = (k_scale * g * tanh(k_scale * depth)) ** 0.5
        Uscale = (
            w_scale * height / 2 * (cosh(k_scale * (depth + height / 2)) / sinh(k_scale * depth))
        )
        scale = Uscale * 8

        # Plot surface elevation
        (eta_line,) = ax1.plot(x, eta, label=model_name, zorder=1000, lw=2)
        eta_color = eta_line.get_color()

        # Plot colocation points
        if hasattr(wave, "eta") and t == 0:
            ax1.plot(wave.x, wave.eta, "kx", ms=2)

        if plot_quiver:
            # Plot the velocity components
            ax1.quiver(
                X,
                Y,
                U,
                V,
                scale_units="height",
                scale=scale,
                color=eta_color,
                width=3.5e-3,
                label=None,
                alpha=0.6,
            )

            # Plot the velocities in the air
            if plot_air:
                X2, Y2 = makeXY(xvec, evec + height / 8, xvec * 0 + ymax, Ny)
                vel = wave.velocity(X2.ravel(), Y2.ravel(), t)
                U2 = vel[:, 0].reshape(X.shape)
                V2 = vel[:, 1].reshape(X.shape)
                ax1.quiver(
                    X2,
                    Y2,
                    U2,
                    V2,
                    scale_units="height",
                    scale=scale,
                    color=eta_color,
                    width=2.0e-3,
                    label=None,
                    alpha=0.6,
                )

        # Plot the velocities in vertical slices
        for i, ax in enumerate(ax2s):
            xi = length / 2 * i / 4
            e = wave.surface_elevation(xi, t)[0]
            y = numpy.linspace(ymin, e, 1000)
            v = wave.velocity(y * 0 + xi, y, t)
            if i < 5:
                ax.plot(v[:, 0], y, c=eta_color)
                vx_max = max(abs(v[:, 0]).max(), vx_max)
            else:
                ax.plot(v[:, 1], y, c=eta_color)
                vz_max = max(abs(v[:, 1]).max(), vz_max)

            if plot_air:
                y = numpy.linspace(e, ymax, 1000)
                v = wave.velocity(y * 0 + xi, y, t)
                if i < 5:
                    ax.plot(v[:, 0], y, ":", c=eta_color)
                    vx_max = max(abs(v[:, 0]).max(), vx_max)
                else:
                    ax.plot(v[:, 1], y, ":", c=eta_color)
                    vz_max = max(abs(v[:, 1]).max(), vz_max)

        # Crest and trough velocities
        eta0 = wave.surface_elevation([0.0, length / 2], t=0)
        u_crest = wave.velocity(0.0, eta0[0], t=0)[0, 0]
        u_trough = wave.velocity(length / 2, eta0[1], t=0)[0, 0]

        # Print some info
        print(
            info_format
            % (
                model_name,
                wave.c,
                2 * pi / (wave.k * wave.c),
                eta.max(),
                eta.min(),
                u_crest,
                u_trough,
                vx_max,
            )
        )

    # Print the velocity scale
    print(head_format % ("scale", "", "", "", "", "", "", "%10.3e" % Uscale))
    if warnings:
        print(warnings)

    # Show a legend if there are more than one wave model being plotted
    if len(model_names) > 1:
        ax1.legend(loc="lower right")
    ax1.set_title("%s waves" % " and ".join(model_names))

    fig2.suptitle("Horisontal (upper) and vertical (lower) velocities")
    ax2s[0].set_title("x = 0 L")
    ax2s[1].set_title("1/8 L")
    ax2s[2].set_title("1/4 L")
    ax2s[3].set_title("3/8 L")
    ax2s[4].set_title("1/2 L")

    # Set axes limits
    ax1.set_ylim(ymin, ymax)
    for i, ax in enumerate(ax2s):
        if i < 5:
            ax.set_xlim(-vx_max, vx_max)
        else:
            ax.set_xlim(-vz_max, vz_max)
        ax.set_xlim(-Uscale * 1.3, Uscale * 1.3)
        if i not in (0, 5):
            ax.set_xticks([])
            ax.set_yticks([])

    fig1.tight_layout()
    fig2.tight_layout(rect=[0, 0.0, 1, 0.95])
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
            Y[i, j] = (1 - F[i, j]) * ymin[j] + F[i, j] * ymax[j]
    return X, Y


def main():
    # Get command line arguments
    import argparse

    parser = argparse.ArgumentParser(prog="raschii.cmd.plot", description="Plot a Raschii wave")
    parser.add_argument(
        "wave_type", help="Name of the wave model. " 'You can specify several, e.g., "Fenton,Airy"'
    )
    parser.add_argument("wave_height", help="Wave height", type=float)
    parser.add_argument("water_depth", help="The still water depth", type=float)
    parser.add_argument("wave_length", help="Distance between peaks", type=float)
    parser.add_argument("-N", type=int, default=10, help="Approximation order")
    parser.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="Allow exceeding breaking criteria",
    )
    parser.add_argument(
        "-v", "--velocities", default=False, action="store_true", help="Show velocities as arrows"
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
        "--blend_height",
        default=-1,
        type=float,
        help="Blend the water and air stream functions a "
        "distance up to improve continuity of velocities",
    )
    parser.add_argument("-t", "--time", default=0.0, type=float, help="The time instance to plot")
    parser.add_argument("--ymin", default=None, type=float, help="Lower vertical axis limit")
    parser.add_argument("--ymax", default=None, type=float, help="Upper vertical axis limit")
    args = parser.parse_args()

    err, warn = check_breaking_criteria(args.wave_height, args.water_depth, args.wave_length)
    if err:
        print(err)
    if warn:
        print(warn)
    if err and not args.force:
        exit(1)

    model_names = args.wave_type.split(",")
    plot_wave(
        model_names=model_names,
        height=args.wave_height,
        depth=args.water_depth,
        length=args.wave_length,
        N=args.N,
        depth_air=args.depth_air,
        blend_height=args.blend_height,
        t=args.time,
        plot_quiver=args.velocities,
        ymin=args.ymin,
        ymax=args.ymax,
    )


if __name__ == "__main__":
    main()
