deltaBS
=======

## A method for quantifying the significance of genetic variation using probabilistic profile-based methods.

The DeltaBS approach uses profile hidden Markov models that take a diverse collection of sequences of the same gene from different organisms and capture patterns of sequence variation that occur commonly in nature. These models are used to score the protein coding genes of different strains of bacteria to produce bitscores, which indicate how well each protein conforms to modelled sequence constraints on that gene. Differences between these bitscores can then be used to identify genes which show a difference in adherence to sequence constraints between strains from different niches. A large difference in bitscore for the same gene may be caused by a change in selective pressures on a gene in a particular niche, leading to increased accumulation of mutations not commonly observed in nature. This accumulation of rare mutations is likely to result in the degradation of the gene, however these mutations could also result in a change to the function or structural stability of the protein. For a given protein family, the difference in bitscore for one strain and the median bitscore for the strain collection (DeltaBS) can be used as an indication of the likelihood that this accumulation of mutations is deleterious.

See the Wiki page for information on setting up and running the software:

https://github.com/UCanCompBio/deltaBS/wiki
