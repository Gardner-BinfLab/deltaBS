deltaBS
=======

Quantifying the significance of genetic variation using probabilistic profile-based methods.

To run, you need to have HMMER3 installed (http://hmmer.janelia.org/). You also need to download the Pfam HMM collection from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release (we use Pfam-A.hmm). When running deltaBS, the option -hd specifies the directory that HMMER3 has been installed in and -hp specifies the directory that the Pfam HMMs have been dowloaded to. 

Once you have installed the Pfam HMMs, hmmpress must be used once to prepare the HMM database for use:

hmmpress [hmmer directory]

To test:

./src/deltaBS.pl -e1 data/tests/esxmtb.embl   -e2 data/tests/esxsmeg.embl  --verbose          -o /tmp -hp /usr/local/bin/ -hd ../delta-bitscore/ -t /tmp

./src/deltaBS.pl -e1 data/tests/FQ312003.embl -e2 data/tests/AE014613.embl --verbose --dirty  -o /tmp -hp /usr/local/bin/ -hd ../delta-bitscore/ -t /tmp

If re-running the analysis, you can skip the hmmscan step by specifying the names of the files that were created in the first run. These follow the format [sequence filename]-pfam_hmmscan1.tbl:
./src/deltaBS.pl -e1 data/tests/FQ312003.embl -e2 data/tests/AE014613.embl --verbose --dirty  -o data/tests/salmonella -hp /usr/local/bin -hd ../delta-bitscore -t data/tests/salmonella -pa1 data/tests/FQ312003.embl-pfam_hmmscan1.tbl -pa2 data/tests/AE014613.embl-pfam_hmmscan1.tbl 


