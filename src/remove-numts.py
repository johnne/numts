#!/usr/bin/env python

import pandas as pd
import sys
from tqdm import tqdm
from argparse import ArgumentParser


def sort_asvs(df):
    return df.sum(axis=1).sort_values(ascending=False)


def get_closest(pids):
    """
    asv: asv under examination
    pids: dataframe of pairwise identities
    """
    if pids.empty:
        return None
    return pids.sort_values("pident", ascending=False).iloc[0].asv2


def main(args):
    if args.asv_list:
        with open(args.asv_list, "r") as fhin:
            subset_asv_list = [x.rstrip() for x in fhin.readlines()]
    countsfile = args.countsfile
    # taxfile = args.taxfile
    idfile = args.idfile
    abundance_factor = args.abundance_factor
    contamination_threshold = args.contamination_threshold
    percent_sample_threshold = args.percent_sample_threshold

    countsdf = pd.read_csv(countsfile, sep="\t", index_col=0)
    sys.stderr.write(f"Counts for {countsdf.shape[0]} ASVs read from {countsfile}")
    # taxdf = pd.read_csv(taxfile, sep="\t", index_col=0)
    # taxdf.loc[taxdf.representative == 1].shape
    pids = pd.read_csv(idfile, sep="\t", header=None, names=["asv1", "asv2", "pident"])
    # remove self-hits
    pids = pids.loc[pids.asv1 != pids.asv2]

    sys.stderr.write("Sorting ASVs by total abundance\n")
    asv_sum = sort_asvs(countsdf)
    asv_list = list(asv_sum.index)
    contamination = {}
    numts = {}

    # get sample where each asv is highest in abundance
    sys.stderr.write("Identifying max-abundance samples for ASVs\n")
    idxmax = countsdf.idxmax(axis=1)

    real_asvs = {
        asv_list[0]: {
            "percent_echo_samples": None,
            "num_echo_samples": None,
            "status": "No1 ASV",
            "n_samples": None,
        }
    }
    if args.asv_list:
        pbar_asv_list = subset_asv_list
    else:
        pbar_asv_list = asv_list[1:]
    pbar = tqdm(pbar_asv_list, desc="checking for numts", unit=" ASVs", ncols=120)
    for asv in pbar:
        pbar.set_postfix(
            {
                "real": len(real_asvs.keys()),
                "numts": len(numts.keys()),
                "contamination": len(contamination.keys()),
            }
        )
        # find the sample in which the ASV is most abundant
        max_sample = idxmax.loc[asv]
        counts = countsdf.loc[:, max_sample]
        asvs_in_sample = list(set(counts.index[counts>0]).difference(numts.keys()))
        _pids = pids.loc[(pids.asv1==asv) & (pids.asv2.isin(asvs_in_sample))]
        # find the closest ASV
        closest_asv = get_closest(_pids)
        # if there are no comparisons for the ASV in this sample (e.g. if no other ASV
        # has a pid above the vsearch minimum pident (84%) then closest_asv is None and we treat this
        # ASV as a real sequence
        if closest_asv is None:
            real_asvs[asv] = {
                "percent_echo_samples": None,
                "num_echo_samples": None,
                "status": "NO PDIST",
                "n_samples": None,
            }
            continue
        # check if closest ASV is more abundant than current ASV
        # if it is not, then the current ASV is not a numt
        if not counts.loc[closest_asv] > counts.loc[asv] * abundance_factor:
            real_asvs[asv] = {
                "percent_echo_samples": None,
                "num_echo_samples": None,
                "status": "MORE ABUNDANT",
                "n_samples": None,
            }
            continue
        ## consistent echo ##
        # get counts for ASV
        sample_counts = countsdf.loc[asv]
        # filter away potential contamination by only checking samples in which the current ASV is
        # above a certain threshold in raw counts
        samples_to_test = sample_counts.loc[
            sample_counts >= contamination_threshold
        ].index
        # if there are no samples in which the current ASV is above the contamination threshold, treat the ASV as a numt
        if len(samples_to_test) == 0:
            contamination[asv] = {
                "nsamples": len(samples_to_test),
                "max_reads": max(sample_counts),
            }
            continue
        # subset to dataframe containing only the ASV and closest ASV, for samples in which the current ASV is present
        comparison_counts = countsdf.loc[[asv, closest_asv], samples_to_test]
        # normalize the dataframe by the counts of the current sequence
        comparison_counts = comparison_counts.div(comparison_counts.loc[asv])
        # now calculate the proportion of samples in which the closest sequence is more abundant (scaled by the abundance factor)
        num_echo_samples = sum(comparison_counts.loc[closest_asv] > abundance_factor)
        percent_echo_samples = num_echo_samples / len(samples_to_test) * 100
        if percent_echo_samples > percent_sample_threshold:
            numts[asv] = {
                "percent_echo_samples": percent_echo_samples,
                "num_echo_samples": num_echo_samples,
                "n_samples": len(samples_to_test),
            }
        else:
            real_asvs[asv] = {
                "percent_echo_samples": percent_echo_samples,
                "num_echo_samples": num_echo_samples,
                "status": "NOT ECHO",
                "n_samples": len(samples_to_test),
            }
    numts_df = pd.DataFrame(numts).T
    numts_df = numts_df.assign(
        status=pd.Series(["NUMT"] * numts_df.shape[0], index=numts_df.index)
    )
    real_df = pd.DataFrame(real_asvs).T
    df = pd.concat([numts_df, real_df])
    df.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser(
        """
    This script identifies NUMTs by taking read counts distribution and pairwise identities into account

    python remove-numts.py counts.tsv pairwise_ids.tsv > numts_results.tsv
    """
    )
    parser.add_argument(
        "countsfile",
        type=str,
        help="Tab-separated file with counts for ASVs (rows) in each sample (columns)",
    )
    parser.add_argument(
        "idfile",
        type=str,
        help="Tab-separated file with pairwise identities for ASVs: <asv1><TAB><asv2><TAB><identity>",
    )
    # parser.add_argument(
    #    "taxfile",
    #    type=str,
    #    help="Tab-separated file with taxonomic assignments for ASVs. Each row should contain one ASV "
    #    "and columns should correspond to relevant taxonomic levels",
    # )
    parser.add_argument(
        "-p",
        "--percent_sample_threshold",
        type=float,
        default=80,
        help="For a sequence under investigation, the closest more abundant sequence has to be more "
        "abundant in this percentage of samples in which both sequences occurr (default: 80)",
    )
    parser.add_argument(
        "-c",
        "--contamination_threshold",
        type=int,
        default=3,
        help="Treat sequences with less than this number of reads (summed across all samples) as contamination."
        " These sequences will be removed prior to numts identification (default: 3)",
    )
    parser.add_argument(
        "-a",
        "--abundance_factor",
        type=float,
        default=1,
        help="In order to be identified as a numt, the closest more abundant sequence has to be more "
        "abundant by this factor (default: 1)",
    )
    parser.add_argument(
        "--asv_list",
        type=str,
        help="(Testing purposes) Supply list of ASVs to check for numts. Will only run the numt checking on these ASVs",
    )

    args = parser.parse_args()
    main(args)
