#!/usr/bin/gawk -f
(/^M/ && $2>1 && $7>-0.01 && $6 > 5 && $10 ~ /[KR][.].+[KR][.]./) {
  hit = spek " " $2 " " $6 " " $7;
  peaks=$8/$9
}
(/^M/ && $2==1 && hit && $7<-0.01 && ($8/$9<peaks)) {
  hit = hit " " $6 " " $7;
  print hit;
}

/^S/ {spek=$2;hit="";}
