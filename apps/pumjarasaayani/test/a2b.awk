# convert from angstroms to Bohrs in a structure file
{ for (i = 1; i <= NF; i++) if ( $i ~ /[\-0-9]/ ) $i = $i * 1.889725989; print }
