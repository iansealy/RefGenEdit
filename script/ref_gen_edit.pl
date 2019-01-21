#!/usr/bin/env perl

# PODNAME: ref_gen_edit.pl
# ABSTRACT: Reference Genome Editor

## Author         : iansealy
## Maintainer     : iansealy
## Created        : 2019-01-14
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;

use Getopt::Long;
use Pod::Usage;
use File::Spec;
use Path::Tiny;
use RefGenEdit qw(edit_genome);

# Default options
my $genome;
my $vcf;
my @sample;
my $method = 'alt';
my $filter;
my $seed;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Ensure reproducible edits if seed set
if ( defined $seed ) {
    srand $seed;
}

edit_genome(
    {
        genome => $genome,
        vcf    => $vcf,
        sample => \@sample,
        method => $method,
        filter => $filter,
    }
);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'genome=s'     => \$genome,
        'vcf=s'        => \$vcf,
        'sample=s@{,}' => \@sample,
        'method=s'     => \$method,
        'filter'       => \$filter,
        'seed=i'       => \$seed,
        'help'         => \$help,
        'man'          => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

__END__

=head1 USAGE

    ref_gen_edit.pl
        [--genome file]
        [--vcf file]
        [--sample sample...]
        [--method alt|other]
        [--filter]
        [--seed seed]
        [--help]
        [--man]

=head1 DESCRIPTION

A script for editing a reference genome using a VCF file.

Currently only SNPs are used when editing. Indels and other variants are
ignored.

Memory usage should be fairly low, provided the genome and VCF are sorted in
the same order.

=head1 EXAMPLES

perl ref_gen_edit.pl --genome genome.fa --vcf snps.vcf

perl ref_gen_edit.pl --genome genome.fa --vcf snps.vcf --sample Sample1 Sample2 --method other --filter --seed 1

=head1 OPTIONS

=over 8

=item B<--genome FILE>

FASTA file containing reference genome.

=item B<--vcf FILE>

VCF file containing variants for editing genome.

=item B<--sample SAMPLE...>

Samples from the VCF file to use when editing genome. By default all samples
are used.

=item B<--method alt|other>

The editing method to use. The default is "alt" and changes the reference to
match one of the ALT alleles in the VCF file. If there are multiple ALT alleles
then one is chosen randomly. The "other" method changes the reference to a base
that doesn't match either the reference allele or one of the ALT alleles. Again,
if there are multiple options then one is chosen randomly.

=item B<--filter>

Only include variants whose FILTER field is set to PASS in the VCF.

=item B<--seed INT>

Random seed (to get a reproducible edited genome).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
