## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package RefGenEdit;
## VERSION
## use critic

# ABSTRACT: Reference Genome Editor

## Author         : iansealy
## Maintainer     : iansealy
## Created        : 2018-12-19
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;

use Readonly;
use List::MoreUtils qw(any pairwise);
use IO::Uncompress::Gunzip qw($GunzipError);
use Exporter qw( import );
our @EXPORT_OK = qw( %METHOD edit_genome get_chr get_variants_by_chr edit_seq );

# Constants
Readonly our $FASTA_OUTPUT_WIDTH         => 80;
Readonly our $VCF_MANDATORY_COLUMN_COUNT => 9;
Readonly our %METHOD                     => (
    ALT         => 1,    # Edit to match ALT allele
    NOT_REF_ALT => 2,    # Edit to match anything other than REF or ALT alleles
);                       # Random if multiple options

=func edit_genome

  Usage       : edit_genome( { genome => $fasta, vcf => $vcf } )
  Purpose     : Edit and output genome
  Returns     : undef
  Parameters  : Hashref {
                    genome => String (the genome file)
                    vcf    => String (the VCF file)
                    sample => String or Arrayref or undef (the required samples)
                    method => String or undef ("alt" or "other")
                    filter => Boolean or undef (only keep PASS positions)
                    fh     => Filehandle or undef (output filehandle)
                }
  Throws      : If genome or VCF arguments are missing or don't exist
                If method is unknown
  Comments    : The "alt" method is the default, where the reference is changed
                to match the ALT allele (or a random ALT allele if there are
                multiple ALT alleles).
                The "other" method changes the reference to a base that's
                different to both the REF and ALT alleles (again, randomly, if
                there are multiple options)

=cut

sub edit_genome {
    my ($arg_ref) = @_;

    confess 'Genome argument missing' if !defined $arg_ref->{genome};
    confess 'VCF argument missing'    if !defined $arg_ref->{vcf};

    my $genome_fh = _open_genome( $arg_ref->{genome} );
    my $vcf_fh    = _open_vcf( $arg_ref->{vcf} );

    my $method =
       !defined $arg_ref->{method}       ? $METHOD{ALT}
      : lc $arg_ref->{method} eq 'alt'   ? $METHOD{ALT}
      : lc $arg_ref->{method} eq 'other' ? $METHOD{NOT_REF_ALT}
      :                                    undef;
    confess sprintf 'Unknown method argument (%s)', $arg_ref->{method}
      if !defined $method;

    while ( my ( $short_name, $long_name, $seq ) = get_chr($genome_fh) ) {
        my $variants =
          get_variants_by_chr( $vcf_fh, $short_name, $arg_ref->{samples},
            $arg_ref->{filter} );
        foreach my $variant ( @{$variants} ) {
            my ( $ref, $alt, $pos ) = @{$variant};
            edit_seq( \$seq, $ref, $alt, $pos, $method );
        }

        #$seq = multiple_edit_seq($seq, $variants, $method);
        _write_chr( $long_name, $seq, $arg_ref->{fh} );
    }

    close $genome_fh;
    close $vcf_fh;

    return;
}

sub _open_genome {
    my ($genome) = @_;

    confess sprintf 'Genome does not exist (%s)', $genome if !-e $genome;
    my $fh;
    if ( $genome =~ m/[.]gz \z/xms ) {
        $fh = IO::Uncompress::Gunzip->new(
            $genome,
            MultiStream => 1,
            Transparent => 0
        ) or confess sprintf 'gunzip failed (%s): %s', $genome, $GunzipError;
    }
    else {
        open $fh, q{<}, $genome;    ## no critic (RequireBriefOpen)
    }

    return $fh;
}

sub _open_vcf {
    my ($vcf) = @_;

    confess sprintf 'VCF does not exist (%s)', $vcf if !-e $vcf;
    my $fh;
    if ( $vcf =~ m/[.]gz \z/xms ) {
        $fh = IO::Uncompress::Gunzip->new(
            $vcf,
            MultiStream => 1,
            Transparent => 0
        ) or confess sprintf 'gunzip failed (%s): %s', $vcf, $GunzipError;
    }
    else {
        open $fh, q{<}, $vcf;    ## no critic (RequireBriefOpen)
    }

    return $fh;
}

sub _write_chr {
    my ( $name, $seq, $fh ) = @_;

    if ( !defined $fh ) {
        $fh = \*STDOUT;
    }

    $seq =~ s/(.{$FASTA_OUTPUT_WIDTH})/$1\n/xmsg;    # Wrap FASTA

    printf {$fh} ">%s\n%s\n", $name, $seq;

    return;
}

=func get_chr

  Usage       : get_chr($fh)
  Purpose     : Get next chromosome from a genome filehandle
  Returns     : Array of chromsome short name, long name and sequence
  Parameters  : Filehandle (the VCF filehandle)
  Throws      : If genome filehandle is missing
  Comments    : None

=cut

{
    my $prev_header_line;

    sub get_chr {
        my ($fh) = @_;

        confess 'Genome filehandle argument missing' if !defined $fh;

        # Get names from header line
        if ( !$prev_header_line ) {
            $prev_header_line = <$fh>;
        }
        chomp $prev_header_line;
        $prev_header_line =~ s/\A >//xms;
        my $long_name  = $prev_header_line;
        my $short_name = $long_name;
        $short_name =~ s/\A \s+//xms;
        $short_name =~ s/\s+ .*//xms;

        my $seq = q{};
        while ( my $line = <$fh> ) {
            if ( $line =~ m/\A >/xms ) {
                $prev_header_line = $line;
                last;
            }
            chomp $line;
            $line =~ s/\s+//xmsg;
            $seq .= $line;
        }
        if ( !$seq ) {

            # EOF
            undef $prev_header_line;
            return;
        }

        return $short_name, $long_name, $seq;
    }
}

=func get_variants_by_chr

  Usage       : get_variants_by_chr($fh, 1)
  Purpose     : Get all REF, ALT & POS for a chromosome from a VCF filehandle
  Returns     : Arrayref of Arrayrefs of REF, ALT and REF position
  Parameters  : Filehandle (the VCF filehandle)
                String (the required chromosome)
                String or Arrayref or undef (the required samples)
                Boolean or undef (whether to only keep PASS positions)
  Throws      : If VCF filehandle or chromosome arguments are missing
                If VCF has no header
                If consecutive variants with same position have different REFs
  Comments    : The filehandle is only read once, so if other chromosomes 
                appear before the required chromosome then their variants are
                stored until they are required, which increases memory usage.

=cut

{
    my $header_seen = 0;
    my @gt_idxs;
    my %gt_data;

    sub get_variants_by_chr {
        my ( $fh, $required_chr, $samples, $only_pass ) = @_;

        confess 'VCF filehandle argument missing' if !defined $fh;
        confess 'Required chromosome argument missing'
          if !defined $required_chr;

        $samples = _convert_to_array($samples);

        while ( my $line = <$fh> ) {

            # Parse header (if encountered)
            if ( $line =~ m/\A[#][#]/xms ) {

                # Ignore meta-information lines, but reinitialise
                $header_seen = 0;
                @gt_idxs     = ();
                %gt_data     = ();
                next;
            }
            if ( $line =~ m/\A[#]CHROM/xms ) {
                $header_seen = 1;
                @gt_idxs     = _parse_header( $line, $samples );
                %gt_data     = ();
                next;
            }
            confess 'VCF has no header' if !$header_seen && @{$samples};

            # Parse and potentially store variant
            my ( $chr, $pos, $ref, $alt ) =
              _parse_line( $line, \@gt_idxs, $only_pass );
            next if !$chr;
            _store_variant( $chr, $pos, $ref, $alt, \%gt_data );

            # Pause parsing variants if already got required chromosome
            last if $chr ne $required_chr && exists $gt_data{$required_chr};
        }

        my $variants = [];
        if ( exists $gt_data{$required_chr} ) {
            $variants = $gt_data{$required_chr};
            delete $gt_data{$required_chr};    # Save memory
        }
        return $variants;
    }
}

sub _parse_header {
    my ( $line, $samples ) = @_;

    chomp $line;
    my @columns = split /\t/xms, $line;

    # Remove mandatory columns and get index for each sample
    # (If no FORMAT or genotypes then indexes will be empty)
    splice @columns, 0, $VCF_MANDATORY_COLUMN_COUNT;
    my @gt_idxs = 0 .. scalar @columns - 1;
    if ( @{$samples} ) {

        # Restrict to just specified samples if necessary
        my %idx_for_sample;
        my $idx = 0;
        foreach my $sample (@columns) {
            $idx_for_sample{$sample} = $idx++;
        }
        @gt_idxs = map { $idx_for_sample{$_} } @{$samples};
        if ( scalar grep { !defined } @gt_idxs ) {
            my @unknown = map { !defined } @gt_idxs;
            my @unknown_samples =
              grep { $_ } pairwise { return $a if $b } @{$samples}, @unknown;
            confess sprintf 'Unknown sample (%s)', join ', ',
              sort @unknown_samples;
        }
    }

    return @gt_idxs;
}

sub _parse_line {
    my ( $line, $gt_idxs, $only_pass ) = @_;

    chomp $line;
    my ( $chr, $pos, undef, $ref, $alt, undef, $filter, undef, undef, @fields )
      = split /\t/xms, $line;

    return () if $only_pass && $filter ne 'PASS';
    return () if length $ref > 1;    # Ignore multiple base reference

    if ( !@{$gt_idxs} ) {

        # Default to checking all genotypes (if present)
        $gt_idxs = [ 0 .. scalar @fields - 1 ];
    }

    # Get required genotypes
    my @gts;
    foreach my $idx ( @{$gt_idxs} ) {
        my $gt = $fields[$idx];
        $gt =~ s/:.*//xms;
        push @gts, $gt;
    }

    # Get ALT SNP alleles
    my %idx_to_alt;
    my $idx = 0;
    foreach my $allele ( split /,/xms, $alt ) {
        $idx++;
        if ( $allele =~ m/\A [AGCT] \z/xms ) {
            $idx_to_alt{$idx} = $allele;
        }
    }
    return () if !%idx_to_alt;    # No ALT SNP alleles

    # Get required ALT SNP alleles
    my %required;
    foreach my $gt (@gts) {
        my @alleles = split /[\/|]/xms, $gt;
        foreach my $required_allele (@alleles) {
            next if $required_allele eq q{0} || $required_allele eq q{.};
            next if !exists $idx_to_alt{$required_allele};
            $required{ $idx_to_alt{$required_allele} } = 1;
        }
    }
    if ( !@gts ) {

        # No FORMAT or genotypes, so all ALT SNP alleles are required
        %required = map { $_ => 1 } values %idx_to_alt;
    }
    return () if !%required;    # No required ALT SNP alleles
    $alt = join q{,}, sort keys %required;

    return $chr, $pos, $ref, $alt;
}

sub _store_variant {
    my ( $chr, $pos, $ref, $alt, $gt_data ) = @_;

    if ( exists $gt_data->{$chr} && $gt_data->{$chr}->[-1]->[-1] == $pos ) {
        confess sprintf 'Different REF at same position (%s:%d)', $chr, $pos
          if $gt_data->{$chr}->[-1]->[0] ne $ref;

        # Merge both ALT
        my %alt = map { $_ => 1 } split /,/xms, $alt;
        foreach my $alt ( split /,/xms, $gt_data->{$chr}->[-1]->[1] ) {
            $alt{$alt} = 1;
        }
        $gt_data->{$chr}->[-1]->[1] = join q{,}, sort keys %alt;
    }
    else {
        push @{ $gt_data->{$chr} }, [ $ref, $alt, $pos ];
    }

    return;
}

=func edit_seq

  Usage       : edit_seq($seq_ref, 'A', 'T', 20, $METHOD{ALT})
  Purpose     : Edit a sequence string at a specific position
  Returns     : undef
  Parameters  : Scalarref - String (the sequence)
                String (the REF base)
                String or Arrayref (the ALT base or bases)
                Int (the REF position)
                Constant ($METHOD{ALT} or $METHOD{NOT_REF_ALT})
  Throws      : If any arguments are missing
                If method is unknown
                If position not a valid integer
                If REF contains multiple bases
                If REF not one of A, G, C, T or N
                If REF does not match genome at specified position
  Comments    : None

=cut

sub edit_seq {
    my ( $seq_ref, $ref, $alt, $pos, $method ) = @_;

    confess 'Sequence argument missing' if !defined $seq_ref;
    confess 'REF argument missing'      if !defined $ref;
    confess 'ALT argument missing'      if !defined $alt;
    confess 'Position argument missing' if !defined $pos;
    confess 'Method argument missing'   if !defined $method;

    confess sprintf 'Unknown method argument (%s)', $method
      if !any { $_ == $method } values %METHOD;

    confess sprintf 'Invalid position (%s)', $pos
      if $pos !~ m/\A \d+ \z/xms || $pos == 0 || $pos >= length ${$seq_ref};

    confess sprintf 'Multiple base REF not yet supported (%s)', $ref
      if length $ref > 1;

    confess sprintf 'REF can only be A, G, C, T or N (%s)', $ref
      if $ref !~ m/[AGCTN]/xmsi;

    my $seq_at_pos = substr ${$seq_ref}, $pos - 1, 1;
    confess sprintf 'REF does not match genome (%s not %s at %d)',
      uc $seq_at_pos, $ref, $pos
      if uc $seq_at_pos ne $ref;

    return if $ref =~ m/N/xmsi;    # No edit if REF is N

    # Convert ALT to array and filter for SNPs
    $alt = _convert_to_array($alt);
    @{$alt} = grep { m/\A [AGCT] \z/xms } @{$alt};    # Ignore indels, etc...

    my @replace_options = @{$alt};                    # Substitute with any ALT
    if ( $method == $METHOD{NOT_REF_ALT} ) {

        # Substitute with any allele that's not REF or ALT
        my %replace_options = map { $_ => 1 } qw(A G C T);
        foreach my $allele ( $ref, @{$alt} ) {
            delete $replace_options{$allele};
        }
        @replace_options = sort keys %replace_options;
    }
    return if !@replace_options;
    my $replacement = $replace_options[ rand scalar @replace_options ];
    if ( $seq_at_pos ne $ref ) {

        # Preserve case of reference
        $replacement = lc $replacement;
    }
    substr ${$seq_ref}, $pos - 1, 1, $replacement;

    return;
}

sub _convert_to_array {
    my ($input) = @_;

    if ( !defined $input ) {
        $input = [];
    }
    if ( ref $input ne 'ARRAY' ) {
        $input = [ split /,/xms, $input ];
    }

    return $input;
}

1;
