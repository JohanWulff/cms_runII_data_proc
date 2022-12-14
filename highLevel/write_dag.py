import json
from glob import glob
import os
from argparse import ArgumentParser
from subprocess import Popen, PIPE


def make_parser():
    parser = ArgumentParser(description="Submit processing of LLR \
Samples")
    parser.add_argument("-s", "--submit_base", type=str, 
                        help="Base dir to submit from")
    parser.add_argument("-o", "--output_dir" ,type=str,
                        help="Dir to write output files to")
    parser.add_argument("-c", "--channel", type=str,
                        help="Channel. Can be either tauTau, muTau, or eTau.")
    parser.add_argument("-j", "--json", type=str,
                        help="JSON File containing paths to samples")
    return parser
                            

def checkmake_dir(path):
    if not os.path.exists(path):
        print(f"{path} does not exist.")
        print("Shall I create it now?")
        yn = input("[y/n] ?")
        if yn.strip().lower() == 'y':
            print('Creating dir(s)!')
            os.makedirs(path)
        else:
            raise ValueError(f"{path} does not exist")


def return_subfile(base_dir, executable):
    arguments = f"-i $(INFILES) -o $(OUTDIR) -e $(EXE) -s $(SAMPLE) \
-y $(YEAR) -c $(CHANNEL)"
    file_str = f"executable={executable}\n\
should_transfer_files = YES\n\
when_to_transfer_output = ON_EXIT\n\
\n\
output                = {base_dir}/out/$(ClusterId).$(ProcId).out\n\
error                 = {base_dir}/err/$(ClusterId).$(ProcId).err\n\
log                   = {base_dir}/log/$(ClusterId).$(ProcId).log\n\
\n\
+JobFlavour = \"longlunch\"\n\
Arguments = {arguments}\n\
queue"
    return file_str


def main(submit_base_dir: str, outdir: str, channel: str, sample_json: str):
    executable = "/eos/user/j/jowulff/res_HH/giles_data_proc/\
CMSSW_10_2_15/src/cms_runII_data_proc/highLevel/executable.py"
    shellscript = "/eos/user/j/jowulff/res_HH/giles_data_proc/\
CMSSW_10_2_15/src/cms_runII_data_proc/highLevel/executable.sh"
    # check if it starts with /afs
    if not submit_base_dir.startswith("/afs"):
        raise ValueError("Submission must happen from /afs!")
    checkmake_dir(submit_base_dir)
    checkmake_dir(outdir)
    # copy executables to /afs. Condor cannot access /eos at the time of writing
    for script in [shellscript, executable]:
        prcs = Popen(f"cp {script} {submit_base_dir}", 
                                shell=True, stdin=PIPE, stdout=PIPE,
                                encoding='utf-8')
        out, err = prcs.communicate()
        if err:
            print(err)
            raise ValueError(f"Unable to move {script} to {submit_base_dir}")

    afs_exe = submit_base_dir+"/executable.py"
    afs_shscript = submit_base_dir+"/executable.sh"

    with open(sample_json) as f:
        d = json.load(f)

    n_files = 0
    for i, sample in enumerate(d):
        print(f"Creating submission dir and writing dag \
files for sample ({i+1}/{len(d)})\r", end="")
        # skip mc samples for now
        if d[sample]["sample_id"] != 0:
            continue
        if not os.path.exists(outdir.rstrip("/")+f"/{sample}"):
            os.mkdir(outdir.rstrip("/")+f"/{sample}")
        submit_dir = submit_base_dir.rstrip("/")+f"/{sample}"
        if not os.path.exists(submit_dir):
            os.mkdir(submit_dir)
            os.mkdir(submit_dir+"/err")
            os.mkdir(submit_dir+"/log")
            os.mkdir(submit_dir+"/out")
        dagfile = submit_dir+f"/{sample}.dag"
        submitfile = submit_dir+f"/{sample}.submit"
        path = d[sample]["path"]
        goodfile = path+"/goodfiles.txt"
        if not os.path.exists(goodfile):
            raise ValueError(f"{sample} does not have a goodfile.txt at \
{path}")
        year = d[sample]["year"]
        with open(goodfile) as gfile:
            gfiles = sorted([line.rstrip() for line in gfile])
            if len(gfiles) == 0:
                print(f"Found {len(gfiles)} files in {goodfile}")
                continue
            n_files+=len(gfiles)
            filechunks = [gfiles[i:i+500] for i in range(0, len(gfiles), 500)]
        with open(dagfile, "x") as dfile:
            for chunk in filechunks:
                print(f"JOB {chunk[0]} {submitfile}", file=dfile)
                print(f'VARS {chunk[0]} INFILES="{" ".join(chunk)}" \
OUTDIR="{outdir.rstrip("/")+f"/{sample}"}" EXE="{afs_shscript}" SAMPLE="{sample}" \
YEAR="{year}" CHANNEL="{channel}"', file=dfile)
        submit_string = return_subfile(base_dir=submit_dir, 
                                       executable=afs_exe)
        with open(submitfile, "x") as subfile:
            print(submit_string, file=subfile)
    print(f"In total {n_files} .root files written to {len(d)} .dag files")


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    main(submit_base_dir=args.submit_base,
         outdir=args.output_dir,
         channel=args.channel,
         sample_json=args.json)