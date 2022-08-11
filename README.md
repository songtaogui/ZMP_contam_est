## Installation

1. Install the listed tools and make sure they could be found in PATH env
2. clone the repos and cd to `/bin` dir
3. make sure the `ZMP_contam_est.sh` is executable, (`chmod +x ZMP_contam_est.sh` would help)
4. run ./ZMP_contam_est.sh to see the usage
5. you may also put the `/path/to/ZMP_contam_est/bin` dir into PATH env for globally executing

## Usage

```bash
------------------------------------------------------------
BlastNT based contamination estimation
------------------------------------------------------------
Dependency: blastn, csvtk, taxonkit, blastdbcmd
------------------------------------------------------------
USAGE:
    bash ZMP_contam_est.sh [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                          show help and exit.
    -t, --threads       <num>   [O]     set threads (default: 2)
    -i, --in            <str>   [R]     input fastq reads, could be gzipped
    -d, --db            <str>   [R]     blast ncbi nt database
    -l, --taxid_list    <str>   [R]     taxid list file, one taxid per line,
        duplicated accessions will be assigned to the first taxid it belongs
        to. Any char following a '#' will be ignored. For example, with file:
            40674    #mammalia
            2759     #eukaryota
        the stats for eukaryota means any eukaryota but the mammalia.
    -n, --num           <int>   [O]     using this number of reads for estimating,
                                        default: 10000
    --evalue            <0-1>   [O]     set e-value threshold for blast (default: 1e-5)
    --identity          <0-1>   [O]     set identity threshold for blast (default: 0.5)
    --acc2taxid         <str>   [O]     set acc2taxid file for blast, format(TSV):
                                        <accession>    <taxid>
                                        default: generated on the fly, which is slow.
    -o, --output        <str>   [O]     output file prefix
    --quiet                             keep quiet, only output fatal errors
    --verbose                           be verbose, output detailed logs
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
```

## Output

- `OUTPUT.subset.fa`: the subset sequences for blast, in plain fasta format
- `OUTPUT.blast`: the blast to nt database results, in tbl format (outfmt 6)
- `OUTPUT.stats.list.gz`: the statistic details for each sequence, with columns of:

    | col | means     | description                                                                                                                                  |
    | --- | --------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
    | 1   | TargetID  | accession ID of the best hit in nt database for the read                                                                                     |
    | 2   | QueryID   | the ID of the query sequence                                                                                                                 |
    | 3   | GroupInfo | the group information, based on the taxids provided through `--taxid_list`, with `:` separated groupOrder, taxid, species name, respectively |
    | 4   | Taxid     | the taxid of the Target sequence                                                                                                             |
- `OUTPUT.stats_species.tsv`: the stats of total read numbers for each target taxid, reverse sorted according to the reads number
- `OUTPUT.stats_group.tsv`: the stats of total read numbers for each taxid listed from `--taxid_list`, reverse sorted according to the reads number

