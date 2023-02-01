from argparse import ArgumentParser
import os
from subprocess import Popen, PIPE
import json
import sys


def make_parser():
    parser = ArgumentParser(description="Resubmit dags that failed on some jobs.\
This script checks if the number of files in the output dir matches the number\
of files in the input_dir.")
    parser.add_argument('-o', '--output_base_dir',
                        type=str, help='output directory.')
    parser.add_argument('-s', '--submit_base_dir', type=str,
                        help='directory to submit from')
    parser.add_argument('-j', '--json', type=str, help='json to read the dasgoclient\
datsets and the names of the subdirs from')
    parser.add_argument("-r", "--resubmit", action='store_true',
                        help="If set, resubmit datasets which aren't complete")
    return parser


def run_dasgoclient(dataset: str):
    """
    Runs dasgoclient and returns a list of files for a given dataset
    """
    cmd = f'dasgoclient -query="file dataset={dataset}"'
    process = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, encoding='utf-8')
    out, err = process.communicate()
    if err:
        print(err)
        sys.exit(1)
    else:
        return out.split()


def n_rootfiles_in_dir(directory: str):
    """
    returns the number of .root files within a given dir
    """
    cmd = f'find {directory} -type f -name "*.root" | wc -l'
    process = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, encoding='utf-8')
    out, err = process.communicate()
    if err:
        print(err)
        sys.exit(1)
    else:
        return int(out)


def submit(dagfile: str):
    if not os.path.exists(dagfile):
        raise ValueError(f"{dagfile} does not exist.")
    else:
        cmd = f"condor_submit_dag {dagfile}"
        print(f"running: {cmd}")
        proc = Popen(cmd, shell=True,
                     stdout=PIPE,
                     stderr=PIPE,
                     encoding='utf-8')
        out, err = proc.communicate()
        if err:
            print(err)
        else:
            print(out)
            print("\n -------------- \n")
            print("success")


def main():
    parser = make_parser()
    args = parser.parse_args()
    output_base_dir = args.output_base_dir
    submit_base_dir = args.submit_base_dir
    resubmit = args.resubmit
    f = open(args.json)
    samples = json.load(f)["2018"]
    for sample in samples:
        n_root = n_rootfiles_in_dir(directory=output_base_dir+f'/{sample}')
        input_path = samples[sample]['Path']
        n_das = n_rootfiles_in_dir(directory=input_path)
        if n_root == n_das:
            print(f"processing complete for {sample}!\r", end="")
        else:
            print(f"number of files in {output_base_dir}/{sample}: {n_root}")
            print(f"number of files in {input_path}: {n_das}")
            if resubmit is True:
                print("Numbers don't match.. Resubmitting")
                dagfile = f"{submit_base_dir}/{sample}/{sample}.dag"
                submit(dagfile=dagfile)
            else:
                print(f"Numbers don't match for {sample}!")


if __name__ == "__main__":
    main()
