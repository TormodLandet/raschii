from raschii import get_wave_model, check_breaking_criteria


def write_swd(swd_file_name, model_name, height, depth, length, N, dt, tmax):
    """
    Write an SWD file for the wave with the given parameters
    """
    WaveClass, _AirClass = get_wave_model(model_name)
    args = dict(height=height, depth=depth, length=length)
    if "N" in WaveClass.required_input:
        args["N"] = N
    wave = WaveClass(**args)
    if wave.warnings:
        print("WARNINGS for %s:\n%s" % (model_name, wave.warnings))

    status = wave.write_swd(swd_file_name, dt=dt, tmax=tmax)
    if status is not NotImplemented:
        print("WRITE SWD DONE\nWrote", swd_file_name)


def main():
    # Get command line arguments
    import argparse

    parser = argparse.ArgumentParser(
        prog="raschii.cmd.swd", description="Write a Raschii wave to file (SWD format)"
    )

    parser.add_argument("swd_file", help="Name of the SWD file to write.")
    parser.add_argument("wave_type", help="Name of the wave model.")
    parser.add_argument("wave_height", help="Wave height", type=float)
    parser.add_argument("water_depth", help="The still water depth", type=float)
    parser.add_argument("wave_length", help="Distance between peaks", type=float)
    parser.add_argument("-N", type=int, default=10, help="Approximation order")
    parser.add_argument("--dt", type=float, default=0.01, help="Timestep")
    parser.add_argument("--tmax", type=float, default=10.0, help="Duration")
    parser.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="Allow exceeding breaking criteria",
    )
    args = parser.parse_args()

    err, warn = check_breaking_criteria(args.wave_height, args.water_depth, args.wave_length)
    if err:
        print(err)
    if warn:
        print(warn)
    if err and not args.force:
        exit(1)

    write_swd(
        swd_file_name=args.swd_file,
        model_name=args.wave_type,
        height=args.wave_height,
        depth=args.water_depth,
        length=args.wave_length,
        N=args.N,
        dt=args.dt,
        tmax=args.tmax,
    )


if __name__ == "__main__":
    main()
