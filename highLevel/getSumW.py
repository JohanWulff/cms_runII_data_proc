import uproot
from argparse import ArgumentParser
import os
import json
from collections import defaultdict

def make_parser():
    parser = ArgumentParser(description="Generate json for submission")
    parser.add_argument("-p", "--paths", type=str, nargs="*",
                        help="Paths containing .root files and goodfiles.txt")
    parser.add_argument("-y", "--year", type=int, help="2016, 17, or 18")
    parser.add_argument("-j", "--json", type=str, required=False, default="",
                        help=".json file to append to. If none, new one will be written.")
    return parser


def main(paths: list, year: int, jsonfile:str = ""):
    if not jsonfile== "":
        with open(jsonfile) as f:
            d = json.load(f)
    else:
        jsonfile = f"{year}.json"
        d = defaultdict(lambda: defaultdict(dict))
    paths = [i for i in paths if os.path.isdir(i)]
    for i, path in enumerate(paths):
        # remove trailing / if present
        path = path.rstrip("/")
        sample = path.split("/")[-1]
        goodfile = path+"/goodfiles.txt"
        if not os.path.exists(goodfile):
            raise ValueError(f"{path} does not have a goodfile.txt")
        n_files = 0
        sum_w = 0
        with open(goodfile) as gfile:
            gfiles = sorted([line.rstrip() for line in gfile])
            if len(gfiles) == 0:
                print(f"Found {len(gfiles)} files in {goodfile}")
                continue
        for j, gfile in enumerate(gfiles):
            print(f"Opening file {j+1}/{len(gfiles)} for sample ({i+1}/{len(paths)}).\r", end="")
            with uproot.open(gfile) as data:
                h_eff = data['h_eff'].to_numpy()
            # h_eff is a tuple like: (array([N1, N2,...], dtype=float32),
            #                         array([0,1,2,3,4,5]))
            sum_w += h_eff[0][0]
        d[year][sample]["Path"] = path
        d[year][sample]["Sum_w"] = sum_w
        n_files+=len(gfiles)
    
    with open(f"{jsonfile}", "w") as f:
        print(f"writing to ./{jsonfile}. You may want to change the name.")
        json.dump(d, f)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    main(paths= args.paths,
         year= args.year,
         jsonfile=args.json)