#!/usr/bin/env bash
# Songtao Gui

# set -o xtrace
# set -o errexit
set -o nounset
set -o pipefail

# >>>>>>>>>>>>>>>>>>>>>>>> Load Common functions >>>>>>>>>>>>>>>>>>>>>>>>
export quiet=FALSE
export verbose=TRUE
export common=$(dirname $0)/../src/common.sh
source $common
if [ $? -ne 0 ];then 
    echo -e "\033[31m\033[7m[ERROR]\033[0m --> Cannot load common functions from easybash lib: $common" >&2
    exit 1;
fi
gst_rcd "Common functions loaded"
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<

usage=$(
cat <<EOF
------------------------------------------------------------
BlastNT based contamination estimation
------------------------------------------------------------
Dependency: blastn, csvtk, taxonkit, blastdbcmd
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

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

EOF
)
if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
export threads=2
export num=10000
export in=""
export db=""
export taxid_list=""
export acc2taxid=""
export evalue=1e-5
export identity=0.5
export output=""
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -t|--threads)
            threads=$2
            shift 2
        ;;
        -i|--in)
            in=$2
            shift 2
        ;;
        -d|--db)
            db=$2
            shift 2
        ;;
        -l|--taxid_list)
            taxid_list=$2
            shift 2
        ;;
        -n|--num)
            num=$2
            shift 2
        ;;
        --evalue)
            evalue=$2
            shift 2
        ;;
        --identity)
            identity=$2
            shift 2
        ;;
        -o|--output)
            output=$2
            shift 2
        ;;
        --acc2taxid)
            acc2taxid=$2
            shift 2
        ;;
        --quiet)
            quiet=TRUE
            shift 1
        ;;
        --verbose)
            verbose=TRUE
            shift 1
        ;;
        *) # unknown flag/switch
            UNKOWN_ARGS+=("$1")
            shift
        ;;
    esac
done
if [ "${#UNKOWN_ARGS[@]}" -gt 0 ];then
    echo "[ERROR] --> Wrong options: \"${UNKOWN_ARGS[@]}\"" >&2
    exit 1
fi
unset UNKOWN_ARGS # restore UNKOWN_ARGS params
# ! Check if required vars are legal
check_var_empty output
check_var_numeric num threads identity
check_files_exists $in ${db}.ndb $taxid_list


# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<

# get acc2taxid file if not provided
if [ ! -s "$acc2taxid" ];then
    gst_log "generate accession to taxid file ..."
    if [ ! -s "${output}.acc2taxid.gz" ];then
        blastdbcmd -db $db -entry all -outfmt '%a %T' | pigz -c > acc2taxid.tsv.gz
        if [ $? -ne 0 ];then gst_err "acc2taxid failed: Non-zero exit";rm -f ${output}.acc2taxid.gz; exit 1;fi
    else
        echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove $PWD/${output}.acc2taxid.gz to force re-run." >&2
    fi
    acc2taxid=${output}.acc2taxid.gz
fi

# get first $num reads
gst_log "get first $num reads ..."
if [ ! -s "${output}.subset.fa" ];then
    seqkit head -j $threads -n $num $in | seqkit fq2fa -j $threads -o ${output}.subset.fa
    if [ $? -ne 0 ];then gst_err "subset reads failed: Non-zero exit";rm -f ${output}.subset.fa; exit 1;fi
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove $PWD/${output}.subset.fa to force re-run." >&2
fi

# run blastn
gst_log "run blastn ..."
if [ ! -s "${output}.blast" ];then
    blastn -query ${output}.subset.fa -db $db -out ${output}.blast -outfmt 6 -num_threads $threads -max_target_seqs 1 -evalue $evalue -perc_identity $identity
    if [ $? -ne 0 ];then gst_err "blast failed: Non-zero exit";rm -f ${output}.blast; exit 1;fi
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove $PWD/${output}.blast to force re-run." >&2
fi

# get taxid accessions
gst_log "get taxid accessions ..."
## USAGE:$0 taxid
fmt_taxid () {
    local taxid=$1
    # -r add rank: species, genus, family, order, class, phylum...
    gst_rcd "dealing with taxid: $taxid ..."
    taxonkit list --ids $taxid -n --indent "" | perl -lane '
        $pre=sprintf("%s:%s",$F[0],join("_", @F[1..$#F])) if $. == 1;
        print "$pre\t$F[0]" unless eof;
    '
}
export -f fmt_taxid
## parallel run
if [ ! -s "${output}.taxid" ];then
    gstcat $taxid_list | perl -ple 's/#.*$//g;s/\s+$//g;' > ${output}.taxid_list.tmp &&\
    parallel -j $threads -k fmt_taxid :::: ${output}.taxid_list.tmp > ${output}.taxid.tmp &&\
    gst_rcd "filter taxid <=> accessions ..." &&\
    gstcat $acc2taxid | perl -F"\t" -lane '
        BEGIN{
            $,="\t";
            #~p get taxid to group info hash
            $inputfile="$ENV{output}.taxid.tmp";
            open(IN,"$inputfile") or die("Cannot open file: $inputfile");
            $gn=1;
            while(<IN>){
                chomp;
                ($group,$taxid)=split(/\t/,$_);
                #~p if duplicated, only assign the first one
                $th{$taxid}=sprintf("G%03d:%s",$gn,$group) if ++$h{$taxid} == 1;
                $gn++ if ++$g{$group} == 1;
            }
            close IN;
            #~p get access in blast hash
            $inputfile="$ENV{output}.blast";
            open(IN,"$inputfile") or die("Cannot open file: $inputfile");
            while(<IN>){
                chomp;
                @blast=split(/\t/,$_);
                $bh{$blast[1]}=1;
            }
            close IN;
        }
        if(exists $th{$F[1]} and exists $bh{$F[0]}){
            print $F[0],$th{$F[1]},$F[1];
        }
    ' > ${output}.taxid
    if [ $? -ne 0 ];then gst_err "get taxid failed: Non-zero exit";rm -f ${output}.taxid; exit 1;fi
    check_files_exists ${output}.taxid
    rm ${output}.taxid.tmp ${output}.taxid_list.tmp
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove $PWD/${output}.taxid to force re-run." >&2
fi


# stat blast out
gst_log "stat blast out ..."
if [ ! -s "${output}.stats_done" ];then
    csvtk cut -tTH -j $threads -f 2,1 ${output}.blast | csvtk uniq -tTH -f 1,2 -o ${output}.stats.tmp &&\
    csvtk join -tTH -j $threads -f 1 ${output}.stats.tmp ${output}.taxid -o ${output}.stats.list.gz &&\
    taxonkit lineage -j $threads -i 4 -L -n ${output}.stats.list.gz |\
        csvtk freq -j $threads -tTH -f 3,5 | csvtk sort -tTH -k 2:nr |\
        csvtk add-header -j $threads -tTH -n Group,Species,Count |\
        csvtk sort -tT -k 3:nr -o ${output}.stats_species.tsv &&\
    csvtk summary -tT -g 1 -f 3:sum -w 0 -j $threads ${output}.stats_species.tsv | csvtk sort -tT -k 2:nr -o ${output}.stats_group.tsv
    if [ $? -ne 0 ];then gst_err "stat blast out failed: Non-zero exit"; exit 1;fi
    echo "done" > ${output}.stats_done
    rm ${output}.stats.tmp
else
    echo -e "\033[35m[WARNING]\033[0m --> Already done, skip running. Remove ${output}.stats_done to force re-run." >&2
fi

gst_log "Done! Outputs:
----------------------------------------------------
$(ls -l ${output}.*)
----------------------------------------------------
stats:
$(ls -l ${output}.stats_*.tsv)
----------------------------------------------------
"
