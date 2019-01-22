# RefGenEdit - Reference Genome Editor
======================================

RefGenEdit is a tool for making edits to a reference genome based on known variation.

Reference bias can occur if a sequencing read is more likely to map successfully to the reference allele or to receive a higher quality score after mapping. An edited reference genome, incorporating known variation, can be used to reduce reference bias.

The genome must be in FASTA format and the variants in VCF format.  Currently only SNPs are used. Indels and other variants in the VCF are ignored.

RefGenEdit can be run in two ways:

  1. The "alt" method is the default and changes the reference to match one of the ALT alleles in the VCF file. If there are multiple ALT alleles then one is chosen randomly.

  2. The "other" method changes the reference to a base that doesn't match either the reference allele or any of the ALT alleles. The idea is that all reads are then equally likely (or unlikely) to map, irrespective of whether they include the reference allele or another allele. Again, if there are multiple possible bases then one is chosen randomly.

You can optionally specify specific samples in the VCF from which to deduce the ALT alleles of interest. Other columns in the VCF are ignored. By default, all the samples are used.

A filter option allows only variants whose FILTER field is set to PASS to be included.

To test RefGenEdit, you can download a reference genome and VCF from Ensembl:

:::
wget http://ftp.ensembl.org/pub/release-95/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-95/variation/vcf/tetraodon_nigroviridis/tetraodon_nigroviridis.vcf.gz
:::

An example command line for using RefGenEdit would then be:

:::
ref_gen_edit.pl --genome Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz --vcf tetraodon_nigroviridis.vcf.gz > Tetraodon_nigroviridis_edited.fa
:::

To see other command line options, do:

:::
ref_gen_edit.pl --help
:::

Or:

:::
ref_gen_edit.pl --man
:::

Memory usage should be fairly low, provided the genome and VCF are sorted in the same order. If they aren't then variants are stored in memory until they are required.

For example, the following test took about 25 minutes (on an early 2015 2.9 GHz Intel Core i5 MacBook Pro) and used over 7 GB of RAM:

:::
wget http://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
ref_gen_edit.pl --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz --vcf common_all_20180418.vcf.gz > Homo_sapiens.GRCh38.dna_sm.primary_assembly.common_all_20180418.fa
:::

But if the VCF was sorted in the same order as the genome then less than 2 GB of RAM were used:

:::
gzip -cd common_all_20180418.vcf.gz | sort -k1,1 -k2,2n | bgzip > common_all_20180418.sorted.vcf.gz
ref_gen_edit.pl --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz --vcf common_all_20180418.sorted.vcf.gz > Homo_sapiens.GRCh38.dna_sm.primary_assembly.common_all_20180418.sorted.fa
:::

## Installation

RefGenEdit is written in Perl and can be installed in the same way as any other Perl module.

For example, installing on an Ubuntu system where you have sudo access might look like this:

:::
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential # Install development tools
sudo apt-get install cpanminus # Install cpanm, which is used for installing Perl dependencies
wget https://github.com/iansealy/RefGenEdit/releases/download/v0.1.0/RefGenEdit-0.1.0.tar.gz
tar zxf RefGenEdit-0.1.0.tar.gz
cd RefGenEdit-0.1.0
perl Makefile.PL
sudo cpanm --installdeps . # Install Perl dependencies
make
make test
sudo make install
:::

To test the installation:

:::
wget http://ftp.ensembl.org/pub/release-95/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-95/variation/vcf/tetraodon_nigroviridis/tetraodon_nigroviridis.vcf.gz
ref_gen_edit.pl --genome Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz --vcf tetraodon_nigroviridis.vcf.gz > Tetraodon_nigroviridis_edited.fa
:::

Alternatively, if you don't have sudo access and want to install in your home directory then installation might look like this:

:::
curl -L https://cpanmin.us | perl - App::cpanminus # Install cpanm, which is used for installing Perl dependencies
~/perl5/bin/cpanm --local-lib=~/perl5 local::lib # Install local::lib, which allows installing Perl modules in a specific directory
eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib) # Activate local::lib
wget https://github.com/iansealy/RefGenEdit/releases/download/v0.1.0/RefGenEdit-0.1.0.tar.gz
tar zxf RefGenEdit-0.1.0.tar.gz
cd RefGenEdit-0.1.0
perl Makefile.PL
~/perl5/bin/cpanm --local-lib=~/perl5 --installdeps . # Install Perl dependencies
make
make test
make install
:::

To test the installation:

:::
wget http://ftp.ensembl.org/pub/release-95/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-95/variation/vcf/tetraodon_nigroviridis/tetraodon_nigroviridis.vcf.gz
~/perl5/bin/ref_gen_edit.pl --genome Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz --vcf tetraodon_nigroviridis.vcf.gz > Tetraodon_nigroviridis_edited.fa
:::
