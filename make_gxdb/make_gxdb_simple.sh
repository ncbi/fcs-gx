#!/usr/bin/env bash

set -euo pipefail
shopt -s expand_aliases
alias awkt="awk -v FS='\t' -v OFS='\t'"
alias sortt="sort -t$'\t'"

inp=$1
out_dir=gxdb_out

# in GX the maximum seq-id length is 39, but CDS seq-ids in *cds_from_genomic.fna.gz are longer, e.g.
# lcl|NW_012132914.1_cds_NP_001297063.1_130327
# So we'll truncate them to just the suffix cds_NP_001297063.1_130327
# NB: also handled in make_gxdb.sh
shorten_cds_seq_ids_rx="s/^>lcl\|\S+_cds_/>cds_/"

# input is a 6-column TSV with following columns: tax-id, species, common-name, BLAST-div, gx-div, path-to-fasta
# e.g.
# 5855    Plasmodium vivax    malaria parasite P. vivax    apicomplexans    prst:alveolates    /am/ftp-genomes/all/GCA/900/093/545/GCA_900093545.1_PvT01/GCA_900093545.1_PvT01_genomic.fna.gz

echo "Validating inputs..."
cat $inp | awkt '(NF != 6 || !int($1) || $5 !~ /^[a-z]{4}:/){ print "Invalid data:", $0; exit(1) }'

mkdir -p $out_dir

# ---------------------------------------------------------------------------
echo "Making db.taxa.tsv"
printf "##[[\"GX taxonomic divisions\",1,1]]\n#tax-id\tspecies\tcommon-name\tBLAST-div\tdiv\n" > $out_dir/db.taxa.tsv
cat $inp | cut -f 1-5 | sort -u                                                               >> $out_dir/db.taxa.tsv


# ---------------------------------------------------------------------------
echo "Making seq_id-tax_id.tsv.tmp"
printf "##[[\"GX seq-id to tax-id mapping\",1,1]]\n" > $out_dir/seq_id-tax_id.tsv.tmp
cat $inp | while read line
do
    fasta_path=$(printf "$line" | cut -f 6)
    zcat -f $fasta_path |
        grep -Po '^>\S+' |  # >seq_id
        perl -pae $shorten_cds_seq_ids_rx | 
        awkt -v line="$line" '1{ print substr($1,2),line }' | # seq-id, tax-id ...
        cut -f 1-2 >> $out_dir/seq_id-tax_id.tsv.tmp
done

# ---------------------------------------------------------------------------
echo "Making db.{gxi,gxs}"
cat $inp | cut -f 6 | xargs zcat -f | pv -Wbrat |
    perl -pae $shorten_cds_seq_ids_rx | 
    gx make-db --seq_id-tax_id=$out_dir/seq_id-tax_id.tsv.tmp --taxa=$out_dir/db.taxa.tsv --output=$out_dir/db.gxi

rm $out_dir/seq_id-tax_id.tsv.tmp

ls -ls $out_dir/db.{gxi,gxs,meta.jsonl,seq_info.tsv,taxa.tsv}
gzip -f $out_dir/db.seq_info.tsv
