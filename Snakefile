import os

def mem_allowed(wildcards, threads):
    return max(threads * 6400, 6400)

configfile: "config.yml"

rule all:
    input:
        filt = expand("results/filter/{dataset}/{filter}_filtered.txt",
                        filter = config["filters"], dataset=config["dataset"])



rule orffinder:
    output:
        "results/orffinder/{dataset}.orfs.faa"
    input:
        "data/{dataset}.fasta"
    log:
        "logs/orffinder/{dataset}.log"
    envmodules:
        "bioinfo-tools",
        "ORFfinder/0.4.3"
    shell:
        """
        ORFfinder -in {input} -ml 30 -g 5 -s 2 -n true -strand plus -out {output} -outfmt 0 > {log} 2>&1
        """

rule longest_orfs:
    output:
        faa="results/orffinder/{dataset}.longest.faa",
        txt="results/orffinder/{dataset}.longest.lengths.txt",
    input:
        rules.orffinder.output
    log:
        "logs/orffinder/{dataset}.longest.log"
    run:
        seqs = {}
        from Bio.SeqIO import parse
        import re
        for record in parse(input[0], "fasta"):
            seq_id = re.sub("ORF\d+_(\w+).+",r"\g<1>",(record.id).split("|")[1])
            l = len(record.seq)
            if seq_id not in seqs.keys():
                seqs[seq_id] = {'seq': str(record.seq), 'id': record.id}
            else:
                if l>len(seqs[seq_id]["seq"]):
                    seqs[seq_id] = {'seq': str(record.seq), 'id': record.id}
        with open(output.faa, 'w') as fhout, open(output.txt, 'w') as fhout_txt:
            for seq_id, d in seqs.items():
                fhout.write(f">{seq_id} {d['id']}\n{d['seq']}\n")
                fhout_txt.write(f"{seq_id}\t{len(d['seq'])}\n")

rule get_hmms:
    output:
        expand("resources/hmm/bold.hmm{suff}", suff=["",".h3f",".h3i",".h3m",".h3p"]),
    log:
        "logs/wget/wget.log"
    shadow: "minimal"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        exec &>{log}
        wget https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline/releases/download/v4.3.0/bold.hmm.{{h3f,h3i,h3m,h3p}}
        wget https://github.com/Hajibabaei-Lab/SCVUC_COI_metabarcode_pipeline/releases/download/v4.3.0/bold.hmm
        mv bold.hmm* {params.outdir}
        """

rule hmmscan:
    output:
        "results/hmmscan/{dataset}.out"
    input:
        hmm=rules.get_hmms.output,
        faa=rules.longest_orfs.output.faa
    log:
        "logs/hmmscan/{dataset}.log"
    envmodules:
        "bioinfo-tools",
        "hmmer/3.3.2"
    threads: 4
    params:
        resource_dir=lambda wildcards, input: os.path.dirname(input.hmm[0])
    shell:
        """
        hmmscan --cpu {threads} --tblout {output[0]} {params.resource_dir}/bold.hmm {input.faa} >/dev/null 2>{log}
        """

def filter_input(wildcards):
    if wildcards.filter == "length":
        return rules.longest_orfs.output.txt
    elif wildcards.filter == "bitscore":
        return rules.hmmscan.output[0]

rule filter:
    """
    The cutoff for short/low outliers :
    25th percentile - (1.5 * interquartile range)

    The cutoff for long outliers :
    75th percentile + (1.5 * interquartile range)
    """
    output:
        filt="results/filter/{dataset}/{filter}_filtered.txt",
        outliers="results/filter/{dataset}/{filter}_outliers.txt"
    input:
        filter_input
    params:
        percentile_upper = 75,
        percentile_lower = 25
    run:
        f = wildcards.filter
        import pandas as pd
        from scipy.stats import iqr
        from numpy import percentile
        pandas_params = {
            'usecols': {'length': [0,1], 'bitscore': [2,8]},
            'sep': {'length': '\t', 'bitscore': '\s+'},
            'dtype': {'length': {'ASV': str, 'length': int}, 'bitscore': {'ASV': str, 'bitscore': float}}
        }
        df = pd.read_csv(input[0], sep=pandas_params['sep'][f], index_col=None, header=None,
            names=["ASV",f], usecols=pandas_params['usecols'][f], comment="#",
            dtype=pandas_params['dtype'][f])
        df.set_index("ASV", inplace=True)
        i = iqr(df[f])
        pl = percentile(df[f], params.percentile_lower)
        pu = percentile(df[f], params.percentile_upper)
        cutoff_lower = pl - (1.5 * i)
        cutoff_upper = pu + (1.5 * i)
        low_outliers = list(df.loc[df[f] < cutoff_lower].index)
        if f == "bitscore":
            high_outliers = []
        else:
            high_outliers = list(df.loc[df[f] > cutoff_upper].index)
        outliers = low_outliers+high_outliers
        df_filt = df.loc[(df[f]<=cutoff_upper)&(df[f]>=cutoff_lower)]
        df_filt.to_csv(output.filt, sep="\t")
        df_outliers = df.loc[outliers]
        df_outliers.to_csv(output.outliers, sep="\t")

rule blastn:
    output:
        out = "results/blastn/{dataset}/{dataset}.blastn.out"
    input:
        fasta = "data/{dataset}.fasta"
    log:
        "logs/{dataset}/blastn.log"
    threads: 10
    params:
        blastn_settings = config["blast"]["settings"],
        blast_dir=config["blast"]["dbdir"],
        out = "$TMPDIR/out"
    envmodules:
        "bioinfo-tols",
        "blast/2.14.1+"
    resources:
        runtime = 60 * 24 * 10
    shell:
        """
        blastn -num_threads {threads} -outfmt 6 -db {params.blast_dir}/nt -query {input.fasta} {params.blastn_settings} -out {params.out} 2 > {log}
        mv {params.out} {output.out}
        """

rule blastdbcmd:
    output:
        "results/blastn/{dataset}/{dataset}.hits.tsv"
    input:
        blastout = rules.blastn.output[0]
    params:
        blast_dir=config["blast"]["dbdir"],
        tmpdir="$TMPDIR"
    envmodules:
        "bioinfo-tols",
        "blast/2.14.1+"
    resources:
        runtime = 60 * 10
    shell:
        """
        cut -f2 {input.blastout} | sort -u > {params.tmpdir}/hits
        blastdbcmd -outfmt "%a\t%l\t%t" -db {params.blast_dir}/nt -entry_batch {params.tmpdir}/hits | sed 's/\\t/\t/g' > {output[0]}
        rm {params.tmpdir}/hits
        """

## VSEARCH ALIGNMENTS ##
rule vsearch_align:
    input:
        fasta="data/{dataset}.fasta",
    output:
        dist="results/vsearch/{dataset}/dist.gz",
    log:
        "logs/vsearch/{dataset}_align.log",
    params:
        dist="$TMPDIR/vsearch/{dataset}/dist",
        fasta="$TMPDIR/vsearch/{dataset}/fasta",
        tmpdir="$TMPDIR/vsearch/{dataset}",
        id=config["vsearch"]["id"],
        iddef=config["vsearch"]["iddef"],
        query_cov=config["vsearch"]["query_cov"],
    threads: config["vsearch"]["threads"]
    conda:
        "envs/vsearch.yml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.18.0"
    resources:
        runtime=60 * 24 * 10,
        mem_mb=mem_allowed,
    shell:
        """
        mkdir -p {params.tmpdir}
        cp {input.fasta} {params.fasta}
        vsearch --usearch_global {params.fasta} --db {params.fasta} --self \
            --userout {params.dist} -userfields query+target+id --maxaccepts 0 --maxrejects 0 \
            --id {params.id} --iddef {params.iddef}  --query_cov {params.query_cov} --threads {threads} > {log} 2>&1
        gzip {params.dist}
        mv {params.dist}.gz {output.dist}
        rm -rf {params.tmpdir} 
        """