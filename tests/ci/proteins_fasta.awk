# Writes a FASTA file for the proteins of a tab-delimited percolator input:
# every protein, in the order of its first appearance, as the concatenation
# of its peptides (without flanking residues and modifications). Used by
# option_matrix.sh for the runs that read a FASTA file (-f <fasta>).
#
# usage: awk -f proteins_fasta.awk input.pin > proteins.fasta
BEGIN { FS = "\t" }
NR == 1 {
  for (i = 1; i <= NF; i++) {
    if ($i == "Peptide") pep = i
  }
  next
}
pep > 0 {
  seq = $pep
  sub(/^[^.]*\./, "", seq)
  sub(/\.[^.]*$/, "", seq)
  gsub(/\[[^]]*\]/, "", seq)
  for (i = pep + 1; i <= NF; i++) {
    if (!($i in prot)) {
      prot[$i] = ""
      order[++n] = $i
    }
    if (!(($i, seq) in seen)) {
      seen[$i, seq] = 1
      prot[$i] = prot[$i] seq
    }
  }
}
END {
  for (i = 1; i <= n; i++) {
    printf(">%s\n%s\n", order[i], prot[order[i]])
  }
}
