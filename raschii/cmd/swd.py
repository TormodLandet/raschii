from raschii import get_wave_model, check_breaking_criteria


def write_swd(swd_file_name, model_name, height, depth, length, N, dt, tmax, amp=1, period=None):
    """
    Write an SWD file for the wave with the given parameters.

    Either *length* or *period* must be provided.  Pass ``length=None`` together
    with a ``period`` value to specify the wave by its period instead.
    """
    WaveClass, _AirClass = get_wave_model(model_name)
    args: dict = dict(height=height, depth=depth)
    if length is not None:
        args["length"] = length
    elif period is not None:
        args["period"] = period
    else:
        from raschii import RaschiiError
        raise RaschiiError("Either length or period must be given")
    if "N" in WaveClass.required_input:
        args["N"] = N
    wave = WaveClass(**args)
    if wave.warnings:
        print("WARNINGS for %s:\n%s" % (model_name, wave.warnings))

    wave.write_swd(swd_file_name, dt=dt, tmax=tmax, amp=amp)
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
    parser.add_argument(
        "wave_length",
        nargs="?",
        default=None,
        type=float,
        help="Distance between wave crests. Mutually exclusive with --period / -T.",
    )
    parser.add_argument(
        "-T",
        "--period",
        type=float,
        default=None,
        help="Wave period in seconds (alternative to the positional wave_length argument).",
    )
    parser.add_argument("-N", type=int, default=10, help="Approximation order")
    parser.add_argument("--dt", type=float, default=0.01, help="Timestep")
    parser.add_argument("--tmax", type=float, default=10.0, help="Duration")
    parser.add_argument(
        "--swd-amp",
        type=int,
        choices=[1, 2, 3],
        default=1,
        help=(
            "SWD amp flag. "
            "1 (default): store potential at z=0 (calm surface); "
            "2: store potential on the wavy free surface; "
            "3: store elevation only (no potential)."
        ),
    )
    parser.add_argument(
        "-f",
        "--force",
        default=False,
        action="store_true",
        help="Allow exceeding breaking criteria",
    )
    args = parser.parse_args()

    if args.wave_length is None and args.period is None:
        parser.error("Either the positional wave_length or --period / -T must be given.")
    if args.wave_length is not None and args.period is not None:
        parser.error("wave_length and --period are mutually exclusive.")

    err, warn = check_breaking_criteria(
        args.wave_height, args.water_depth, length=args.wave_length, period=args.period
    )
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
        period=args.period,
        N=args.N,
        dt=args.dt,
        tmax=args.tmax,
        amp=args.swd_amp,
    )


if __name__ == "__main__":
    main()
