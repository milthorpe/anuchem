# convert from angstroms to Bohrs in a structure file
# skips first three fields which are number of atoms, structure name and basis set
{ for (i = 2; i <= NF; i++) if ( $i ~ /^[-]?[0-9]*\.[0-9]*$/ ) $i = $i * 1.889725989; print }
