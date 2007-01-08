#!/usr/bin/gawk -f
(/^M/ && oldPep && ($7>-0.00001) && $10 ~ /[KR][.].+[KR][.]./ && oldPep) {
  id2=$10
  gsub("[IL]","",id2)
  if (id != id2)
    print spek " " $2 " " $6  " " $7 " " $10 " " oldPep
}
(/^M/ && ($7>-0.00001) && $10 ~ /[KR][.].+[KR][.]./) {
  oldPep = $10;
  id=$10;
  gsub("[IL]","",id)
}
/^S/ {spek=$2;oldPep=""}
