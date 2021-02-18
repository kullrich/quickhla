# quickhla

## About

`quickhla` is inspired by [https://github.com/ExpressionAnalysis/HLAProfiler](https://github.com/ExpressionAnalysis/HLAProfiler)
and [https://sourceforge.net/projects/targt-pipeline/](https://sourceforge.net/projects/targt-pipeline/)
to use k-mer content of next generation sequencing reads and subsequent mapping to classify HLA types in a sample.

`quickhla` heavily relies on the data from the [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) database.

![quickhla workflow](quickhla.png)

## Usage

First, the HLA k-mer databases need to be build (can take up to 30 minutes):

```bash
# download hla fasta files from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/
wget ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_nuc.fasta
# build databases with k-mer length of 35
quickhla build -gen hla_gen.fasta -nuc hla_nuc.fasta
```

see help

```
quickhla build -h

usage: quickhla <sub-script> [options] [<arguments>...] build
       [-h] [-gen GEN] [-nuc NUC] [-o O] [-kb KB] [-bb BB] [-hb HB] [-kl KL]
       [-ml ML] [-ms MS]

optional arguments:
  -h, --help  show this help message and exit
  -gen GEN    input gen file [mandatory]
  -nuc NUC    input nuc file [mandatory]
  -o O        output kraken2 database directory [default: hla.db]
  -kb KB      specify kraken2-build binary [if not given assumes to be in
              PATH]
  -bb BB      specify bowtie2-build binary [if not given assumes to be in
              PATH]
  -hb HB      specify hisat2-build binary [if not given assumes to be in PATH]
  -kl KL      kraken2 kmer-len array [default: 35] can be an array of values
              and needs to be same length as -ml and -ms
  -ml ML      kraken2 minimizer-len array [default: 31] if array needs to be
              same length as -kl and -ms
  -ms MS      kraken2 minimizer-space array [default: 7] if array needs to be
              same length as -kl and -ml
```

For paired reads `FORWARD.fq` and `REVERSE.fq` using the standard `nuc` database created for four-digits `4d` , the basic usage is:

```bash
quikhla classify -f FORWARD.fq -r REVERSE.fq -d hla.db -db hla.nuc.4d.35 -n 2
```

see help

```
quickhla classify -h

usage: quickhla <sub-script> [options] [<arguments>...] classify
       [-h] [-f F] [-r R] [-d D] [-db DB] [-o O] [-t T] [-kb KB] [-alg ALG]
       [-bb BB] [-bo BO] [-hb HB] [-ho HO] [-n N]

optional arguments:
  -h, --help  show this help message and exit
  -f F        specify forward fastq [mandatory]
  -r R        specify reverse fastq [mandatory]
  -d D        specify db directory [mandatory]
  -db DB      specify db for classification [mandatory]
  -o O        specify output prefix [default: out]
  -t T        specify number threads [default: 1]
  -kb KB      specify kraken2 binary [if not given assumes to be in PATH]
  -alg ALG    specify aligner [default: hisat2]
  -bb BB      specify bowtie2 binary [if not given assumes to be in PATH]
  -bo BO      specify bowtie2 options [default: --very-fast --no-unal
              --ignore-quals -k 20]
  -hb HB      specify hisat2 binary [if not given assumes to be in PATH]
  -ho HO      specify hisat2 options [default: --fast -k 20]
  -n N        specify number of top hits to report [default: show all]
```

## Output

quickhla outputs two plain text files with the HLA classifications (1. from k-mer and 2. from mapping)

The output are:

## Third-party Dependencies

quickhla requires:

- kraken2 [https://ccb.jhu.edu/software/kraken2/](https://ccb.jhu.edu/software/kraken2/)
- bowtie2 [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- hisat2 [https://daehwankimlab.github.io/hisat2/](https://daehwankimlab.github.io/hisat2/)
- biopython [https://biopython.org/](https://biopython.org/)
- numpy [https://numpy.org/](https://numpy.org/)

You can install them manually (see individual download pages):

- kraken2 [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)
- bowtie2 [https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2)
- hisat2 [https://daehwankimlab.github.io/hisat2/download/](https://daehwankimlab.github.io/hisat2/download/)
- biopython [https://biopython.org/wiki/Download](https://biopython.org/wiki/Download)
- numpy [https://numpy.org/install/](https://numpy.org/install/)

Or you can install them via conda:

- kraken2 [https://anaconda.org/bioconda/kraken2](https://anaconda.org/bioconda/kraken2)
- bowtie2 [https://anaconda.org/bioconda/bowtie2](https://anaconda.org/bioconda/bowtie2)
- hista2 [https://anaconda.org/biobuilds/hisat2](https://anaconda.org/biobuilds/hisat2)
- biopython [https://anaconda.org/anaconda/biopython](https://anaconda.org/anaconda/biopython)
- numpy [https://anaconda.org/anaconda/numpy](https://anaconda.org/anaconda/numpy)

```
conda install -c bioconda kraken2
conda install -c bioconda bowtie2
conda install -c biobuilds hisat2
conda install -c anaconda biopython
conda install -c anaconda numpy
```

## Installation

### via Conda

### PyPI

### Manually

```bash
git clone https://gitlab.gwdg.de/kristian.ullrich/quickhla.git
cd quickhla
python setup.py install
```

### Manually in conda env

```bash
conda create -n quickhla python=3.8
conda activate quickhla
conda install -c bioconda kraken2
conda install -c bioconda bowtie2
conda install -c biobuilds hisat2
conda install -c anaconda biopython
conda install -c anaconda numpy
git clone https://gitlab.gwdg.de/kristian.ullrich/quickhla.git
cd quickhla
python setup.py install
```

## Example

Test HLA typing for individuals from [https://www.internationalgenome.org/](https://www.internationalgenome.org/)

```
# download test data
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/SRR099966_1.filt.fastq.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/sequence_read/SRR099966_2.filt.fastq.gz
# classify reads
quickhla classify -f SRR099966_1.filt.fastq.gz -r SRR099966_2.filt.fastq.gz -d hla.db -db hla.nuc.4d.35 -n 2
```

## Known Issues

## Contributing Code

If you would like to contribute to CRBHits, please file an issue so that one can establish a statement of need, avoid redundant work, and track progress on your contribution.

Before you do a pull request, you should always file an issue and make sure that someone from the CRBHits developer team agrees that it’s a problem, and is happy with your basic proposal for fixing it.

Once an issue has been filed and we've identified how to best orient your contribution with package development as a whole, [fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [main repo](https://gitlab.gwdg.de/kristian.ullrich/quickhla), branch off a [feature branch](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/about-branches) from `master`, [commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/committing-and-reviewing-changes-to-your-project) and [push](https://docs.github.com/en/github/using-git/pushing-commits-to-a-remote-repository) your changes to your fork and submit a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests) for `quickhla:master`.

By contributing to this project, you agree to abide by the Code of Conduct terms.

## Bug reports

Please report any errors or requests regarding [quickhla](https://gitlab.gwdg.de/kristian.ullrich/quickhla) to Kristian Ullrich (ullrich@evolbio.mpg.de)

or use the issue tracker at [https://gitlab.gwdg.de/kristian.ullrich/quickhla/issues](https://gitlab.gwdg.de/kristian.ullrich/quickhla/issues)

## Code of Conduct - Participation guidelines

This repository adhere to [Contributor Covenant](http://contributor-covenant.org) code of conduct for in any interactions you have within this project. (see [Code of Conduct](https://gitlab.gwdg.de/kristian.ullrich/quickhla/-/blob/master/CODE_OF_CONDUCT.md))

See also the policy against sexualized discrimination, harassment and violence for the Max Planck Society [Code-of-Conduct](https://www.mpg.de/11961177/code-of-conduct-en.pdf).

By contributing to this project, you agree to abide by its terms.

## References

1. Buchkovich M.L., Brown C.C., Robasky K., Chai S., Westfall S., Vincent B.G., Weimer E.T., Powers J.G. (2017)
   **HLAProfiler utilizes k-mer profiles to improve HLA calling accuracy for rare and common alleles in RNA-seq data.**
  *Genome Med*, **9, 86**. 
   [https://doi.org/10.1186/s13073-017-0473-6](https://doi.org/10.1186/s13073-017-0473-6)

2. Pierini F., Nutsua M., Böhme L., Özer O., Bonczarowska J., Susat J., Franke A., Nebel A., Krause-Kyora B., Lenz T.L. (2020)
   **Targeted analysis of polymorphic loci from low-coverage shotgun sequence data allows accurate genotyping of HLA genes in historical human populations.**
   *Sci Rep*, **10, 7339**. [https://doi.org/10.1038/s41598-020-64312-w](https://doi.org/10.1038/s41598-020-64312-w)

3. Robinson J., Barker D.J., Georgiou X., Cooper M.A., Flicek P., Marsh S.G.E. (2020)
   **IPD-IMGT/HLA Database**
   *Nucleic Acids Research*, **48(D1), D948–D955**. [https://doi.org/10.1093/nar/gkz950](https://doi.org/10.1093/nar/gkz950)

