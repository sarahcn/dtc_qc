This is an exploration of quality control aspects for potentially clinically relevant SNPs in DTC-GT data.
The impetus was a [finding from Bastian](see (https://twitter.com/gedankenstuecke/status/981342202600161280)
that an APOL1 SNP in a [news story](https://www.huffingtonpost.com/entry/home-genetic-test-false-positives_us_5ac27188e4b04646b6451c42) 
had 13% missingness in openSNP users.

Question: do potentially clinically relevant SNPs in DTC data have a higher missingness compared to other SNPs?
Follow-up question: if yes, could this be due to different MAF profiles (i.e. SNPs with phentoypic impact may be rarer and thus harder to genotype) OR to an intentional data processing step (i.e. masking/setting to missing) on the part of DTC companies?
These are speculations, but as a first pass way to look at this, we:

1. Pulled down the Illumina annotations for Omni Express (a base array used for several years by 23andMe)
1. Used this annotation to identify all missense SNPs in any of the ACMG 59 genes (using missense as an, albeit limited, proxy for potential clinical relevance)
	* for above steps, see Sarah's _dtc_qc_exploration.R_
	* see product in _omniex_missense_acmg.txt_
	
With openSNP data, Bastian can now examine rates of missingess at these SNPs, and also potentially compare annotated MAF (from Illumina) against openSNP users.

SN, 4/28/18