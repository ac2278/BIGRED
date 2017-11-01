# BIGRED

[Ariel W Chan](https://plbrgen.cals.cornell.edu/people/ariel-chan)

---

This file is part of BIGRED.

BIGRED is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BIGRED is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BIGRED.  If not, see <http://www.gnu.org/licenses/>.

---

Researchers typically sequence a given individual multiple times, either re-sequencing the same DNA sample (technical replication) or sequencing different DNA samples collected on the same individual (biological replication) or both. Before merging the data from these replicate sequence runs, it is important to verify that no errors, such as DNA contamination or mix-ups, occurred during the data collection pipeline. BTRED is a R package that performs said task. Using Bayes Theorem, BTRED calculates the posterior probability distribution over the set of relations describing the putative replicates of an individual and infers which of the samples originated from an identical genotypic source. BTRED addresses key limitations of existing methods and produced highly accurate results in simulation experiments. BTRED requires no genotype calling, haplotype phasing, or genotype imputation, making it a suitable tool for studies relying on low- or medium-depth, high-throughput sequence data. 

---

##### Required files:
1. A header file for each chromosome
2. An alias text file for each proband
3. An allelic depth (AD) data file for each proband and chromosome

---
#### 1. A header file for each chromosome (mandatory)
  
**Description:**

Each line in this file corresponds to a site harboring a single nucleotide polymorphism (SNP). 
The file consists of 3 fix, mandatory (tab-delimited) columns. These columns are as follows:

1. REF_FREQ: the estimated frequency of allele A (reference allele) at a given site
2. CHROM: the chromosome number associated with a given site
3. POS: the physical position (in bp) of a given site

**Header Line:**

The file must include a header line, naming the 3 fixed, mandatory columns as ‘REF_FREQ’, ‘CHROM’, and ‘POS’ (no quotes). 

**Missing Values:**

Missing values are specified with ‘NA’ (no quotes). 

**File Naming Format:**

chr[chromosome]_[headersuffix], where [chromosome] must consist of 3 integers. The user should pad [chromosome] with leading zeros if necessary.

**Example:**

We show the header line and the first three rows of file chr001_gbsheader below:

| REF_FREQ   | CHROM     | POS    |
| ---------- | --------- | ------ |
| NA         | 1         | 21594  |
| NA         | 1         | 21597  |
| 0.962656   | 1         | 26576  |

#### 2. An alias text file for each proband (mandatory).
  
**Description:**

This text file lists the sample IDs of the _k_ putative replicates of the proband in question (one sample ID per line). 

**Header Line:**

The file should contain no header line.

**Missing Values:**

The file should contain no missing values. 

**File Naming Format:**

[proband]_aliases.txt

**Example:**

We collected, extracted, and sequenced DNA from proband I011206 a total of _k_ = 3 times, assigning the names 'I011206:250300325', 'I011206:250099287', and 'TMS011206:250134484' to the three samples. We show the contents of file I011206_aliases.txt below:

| I011206:250300325
| I011206:250099287
| TMS011206:250134484

#### 3. An allelic depth (AD) data file containing AD data for the _k_ putative replicates of the proband in question (one file per chromosome)

**Description:**

AD data consists of two (comma-separated) integers, representing the observed counts of allele A (reference allele) and B (alternative allele) at a given site for a given sample. Each line in the file corresponds to a site harboring a SNP. The first 2 columns are fix, mandatory (tab-delimited) columns. These columns are as follows:

1. CHROM
2. POS

Columns 3 through _k_ + 2 list the AD data for the _k_ putative replicates (the IDs of which are listed in the alias file) of the proband in question. 

**Header Line:**

The file must include a header line, naming the first two columns as ‘CHROM’ and ‘POS’ (no quotes) and naming trailing columns of AD data by their corresponding sequence IDs (i.e. the names listed in the alias file).

**Missing Values:**

The file should contain no missing values. If a given sample has been sequenced zero times at a given site, represent this scenario with '0,0' (no quotes) rather than 'NA'.

**File Naming Format:**

chr[chromosome]_[proband].AD.FORMAT, where [chromosome] must consist of 3 integers. The user should pad [chromosome] with leading zeros if necessary.

**Example:**

We show the header line and the first three rows of file chr001_I011206.AD.FORMAT below:

| CHROM | POS   | I011206:250300325 | I011206:250099287 | TMS011206:250134484 |
| ----- | ---   | ----------------- | ----------------- | ------------------- |
| 1	    | 21594	| 2,0	              | 1,0	              | 3,0                 |
| 1	    | 21597	| 2,0	              | 1,0	              | 3,0                 |
| 1	    | 26576	| 5,0	              | 4,0	              | 7,0                 |

Notice how column names for columns 3 through 5 are listed in I011206_aliases.txt. These two files must match.


