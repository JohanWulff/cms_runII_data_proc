#!/usr/bin/bash

usage() { echo "Usage: $0 [-i <input file> ] [-o <outfile>] [-s sample] [-y year] [-c channel]" 1>&2; exit 1; }
while getopts "i:o:s:y:c:" opt; do
    case "$opt" in
        i) INFILE=$OPTARG
            ;;
        o) OUTFILE=$OPTARG
            ;;
        s) SAMPLE=$OPTARG
            ;;
        y) YEAR=$OPTARG
            ;;
        c) CHANNEL=$OPTARG
            ;;
        *)
            echo "Invalid argument $OPTARG";
            usage
            exit 1
    esac
done

CMSSW_SRC="/eos/user/j/jowulff/res_HH/giles_data_proc/CMSSW_10_2_15/src"
EXE="/eos/user/j/jowulff/res_HH/giles_data_proc/CMSSW_10_2_15/bin/slc7_amd64_gcc700/RunLoop"
cd $CMSSW_SRC || exit 1
cmsenv
${EXE} -i $INFILE -o $OUTFILE --sample $SAMPLE -y $YEAR -c $CHANNEL
exit 0