# convert from angstroms to Bohrs in a structure file
awk '{ for (i = 1; i <= NF; i++) if ($i < 0) $i = $i * 1.889725989; print }'
