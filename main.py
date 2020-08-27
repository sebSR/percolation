import simulation.burning as bd
import simulation.hoshen_kopelman as hk
import argparse
import json
import sys


def parameters(filename):
    with open(filename, "r+") as config_file:
        config = json.loads(config_file.read())
    return config


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--method", type=str, help="burning or clustering", required=True)
    args = parser.parse_args()

    filename = "config.json"
    kwargs = parameters(filename)

    if args.method == "burning":
        bd.main_burning(**kwargs)
    elif args.method == "clustering":
        hk.main_hoshen_kopelman(**kwargs)
    else:
        sys.exit("Wrong method")
