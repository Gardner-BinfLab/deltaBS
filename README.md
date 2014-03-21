deltaBS
=======

Quantifying the significance of genetic variation using probabilistic profile-based methods.

Write proper test code, in the mean time:
./src/deltaBS.pl -e1 data/tests/esxmtb.embl   -e2 data/tests/esxsmeg.embl  --verbose          -o /tmp -hp /usr/local/bin/ -hd ../delta-bitscore/ -t /tmp

./src/deltaBS.pl -e1 data/tests/FQ312003.embl -e2 data/tests/AE014613.embl --verbose --dirty  -o /tmp -hp /usr/local/bin/ -hd ../delta-bitscore/ -t /tmp

./src/deltaBS.pl -e1 data/tests/FQ312003.embl -e2 data/tests/AE014613.embl --verbose --dirty  -o data/tests/salmonella -hp /usr/local/bin -hd ../delta-bitscore -t data/tests/salmonella -pa1 data/tests/FQ312003.embl-pfam_hmmscan1.tbl -pa2 data/tests/AE014613.embl-pfam_hmmscan1.tbl 


