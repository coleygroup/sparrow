from sparrow.cli import run, args


def main():
    parser = args
    args = parser.parse_args()
    func = args.func
    del args.func

    func(args)


if __name__ == "__main__":
    main()