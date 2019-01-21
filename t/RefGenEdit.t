use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Capture::Tiny qw(capture_stdout);

use RefGenEdit qw( %METHOD edit_genome get_chr get_variants_by_chr edit_seq );

my @edit_genome_like_data = (
    [
        't/data/test.fa', 't/data/test.vcf', undef,
        ">1 chr1\nTGCTCGCTAGCT\n>2 chr2\nACGATCGATCGACCGA\n",
        'Genome default'
    ],
    [
        't/data/test.fa', 't/data/test.vcf', 'alt',
        ">1 chr1\nTGCTCGCTAGCT\n>2 chr2\nACGATCGATCGACCGA\n",
        'Genome ALT'
    ],
    [
        't/data/test.fa', 't/data/test.vcf', 'other',
        ">1 chr1\n[CG]GCT[GT]GCTAGCT\n>2 chr2\n[CG]CGATCGATCGA[AG]CGA\n",
        'Genome NOT_REF_ALT'
    ],
    [
        't/data/test.fa.gz', 't/data/test.vcf.gz', undef,
        ">1 chr1\nTGCTCGCTAGCT\n>2 chr2\nACGATCGATCGACCGA\n",
        'Compressed genome & VCF'
    ],
    [
        't/data/long.fa',
        't/data/test.vcf',
        undef,
">1 chr1\nTGCTCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\nAGCTAGCTAGCTAGCT\n>2 chr2\nACGATCGATCGACCGA\n",
        '>80 bp chr'
    ],
);

my @edit_genome_throws_data = (
    [ undef,            't/data/test.vcf', undef, 'Genome argument missing' ],
    [ 't/data/test.fa', undef,             undef, 'VCF argument missing' ],
    [ 't/data/missing.fa', 't/data/test.vcf', undef, 'Genome does not exist' ],
    [ 't/data/test.fa',     't/data/missing.vcf', undef, 'VCF does not exist' ],
    [ 't/data/notgz.fa.gz', 't/data/test.vcf',    undef, 'gunzip failed' ],
    [ 't/data/test.fa', 't/data/notgz.vcf.gz', undef, 'gunzip failed' ],
    [
        't/data/test.fa', 't/data/test.vcf',
        'unknown',        'Unknown method argument'
    ],
);

my $genome_string = <<'EOF';
>1 Chr1
AGCT
AGCT
>2 Chr2
TCGA
TCGA
EOF

my @get_chr_is_data = (
    [ $genome_string, q{1},  '1 Chr1', 'AGCTAGCT', 'Genome 1' ],
    [ $genome_string, q{2},  '2 Chr2', 'TCGATCGA', 'Genome 2' ],
    [ $genome_string, undef, undef,    undef,      'Genome EOF' ],
);

my @get_chr_throws_data = ( [ undef, 'Genome filehandle argument missing' ], );

my $vcf_string = <<'EOF';
##fileformat=VCFv4.3
#CHROM POS ID REF ALT  QUAL FILTER INFO FORMAT Sample1 Sample2
1      12  .  A   T    10   PASS   .    GT     0/1     0/0
1      14  .  G   T,C  10   PASS   .    GT     1/1     1/2
1      16  .  A   T    10   PASS   .    GT     0/1     0/0
1      16  .  A   C    10   PASS   .    GT     ./.     0/1
1      18  .  A   TAT  10   PASS   .    GT     0/0     1/1
1      20  .  A   T    1    q10    .    GT     0/0     0/1
1      22  .  AT  TA   10   PASS   .    GT     0/0     1/1
1      24  .  A   T,TA 10   PASS   .    GT     0/2     2/2
2      12  .  A   T    10   PASS   .    GT     0/1     .
EOF

my $vcf_string_noformat = <<'EOF';
##fileformat=VCFv4.3
#CHROM POS ID REF ALT  QUAL FILTER INFO
1      12  .  A   T    10   PASS   .
1      14  .  G   T,C  10   PASS   .
1      16  .  A   T    10   PASS   .
1      16  .  A   C    10   PASS   .
1      18  .  A   TAT  10   PASS   .
1      20  .  A   T    1    q10    .
1      22  .  AT  TA   10   PASS   .
1      24  .  A   T,TA 10   PASS   .
2      12  .  A   T    10   PASS   .
EOF

my $vcf_string_noheader = <<"EOF";
##fileformat=VCFv4.3
1      12  .  A   T   10    PASS   .    GT     0/1     0/0
EOF

my $vcf_string_diff_ref = <<"EOF";
##fileformat=VCFv4.3
1      16  .  A   T   10    PASS   .    GT     0/1     0/0
1      16  .  G   C   10    PASS   .    GT     0/0     0/1
EOF

my @get_variants_by_chr_is_data = (
    [ $vcf_string, q{1}, undef,             undef, 4, 'Chr1' ],
    [ $vcf_string, q{2}, undef,             undef, 1, 'Chr2' ],
    [ $vcf_string, q{3}, undef,             undef, 0, 'Chr3' ],
    [ $vcf_string, q{1}, 'Sample1',         undef, 3, 'Chr1 Sample1' ],
    [ $vcf_string, q{2}, 'Sample1',         undef, 1, 'Chr2 Sample1' ],
    [ $vcf_string, q{3}, 'Sample1',         undef, 0, 'Chr3 Sample1' ],
    [ $vcf_string, q{1}, 'Sample2',         undef, 3, 'Chr1 Sample2' ],
    [ $vcf_string, q{2}, 'Sample2',         undef, 0, 'Chr2 Sample2' ],
    [ $vcf_string, q{3}, 'Sample2',         undef, 0, 'Chr3 Sample2' ],
    [ $vcf_string, q{1}, undef,             1,     3, 'Chr1 PASS' ],
    [ $vcf_string, q{2}, undef,             1,     1, 'Chr2 PASS' ],
    [ $vcf_string, q{3}, undef,             1,     0, 'Chr3 PASS' ],
    [ $vcf_string, q{1}, 'Sample1',         1,     3, 'Chr1 Sample1 PASS' ],
    [ $vcf_string, q{2}, 'Sample1',         1,     1, 'Chr2 Sample1 PASS' ],
    [ $vcf_string, q{3}, 'Sample1',         1,     0, 'Chr3 Sample1 PASS' ],
    [ $vcf_string, q{1}, 'Sample2',         1,     2, 'Chr1 Sample2 PASS' ],
    [ $vcf_string, q{2}, 'Sample2',         1,     0, 'Chr2 Sample2 PASS' ],
    [ $vcf_string, q{3}, 'Sample2',         1,     0, 'Chr3 Sample2 PASS' ],
    [ $vcf_string, q{1}, 'Sample1,Sample2', undef, 4, 'Chr1 Sample1+Sample2' ],
    [ $vcf_string, q{2}, 'Sample1,Sample2', undef, 1, 'Chr2 Sample1+Sample2' ],
    [ $vcf_string, q{3}, 'Sample1,Sample2', undef, 0, 'Chr3 Sample1+Sample2' ],
    [
        $vcf_string, q{1}, [ 'Sample1', 'Sample2' ],
        undef, 4, 'Chr1 [Sample1,Sample2]'
    ],
    [
        $vcf_string, q{2}, [ 'Sample1', 'Sample2' ],
        undef, 1, 'Chr2 [Sample1,Sample2]'
    ],
    [
        $vcf_string, q{3}, [ 'Sample1', 'Sample2' ],
        undef, 0, 'Chr3 [Sample1,Sample2]'
    ],
    [ $vcf_string_noformat, q{1}, undef, undef, 5, 'Chr1 no format' ],
    [ $vcf_string_noformat, q{2}, undef, undef, 1, 'Chr2 no format' ],
    [ $vcf_string_noformat, q{3}, undef, undef, 0, 'Chr3 no format' ],
    [ $vcf_string_noheader, q{1}, undef, undef, 1, 'Chr1 no header' ],
);

no warnings 'qw';
my @get_variants_by_chr_is_deeply_data = (
    [ $vcf_string, q{1}, undef, undef, 0, [qw(A T 12)],   'Chr1:12' ],
    [ $vcf_string, q{1}, undef, undef, 1, [qw(G C,T 14)], 'Chr1:14' ],
    [ $vcf_string, q{1}, undef, undef, 2, [qw(A C,T 16)], 'Chr1:16' ],
    [ $vcf_string, q{1}, undef, undef, 3, [qw(A T 20)],   'Chr1:20' ],
    [ $vcf_string, q{2}, undef, undef, 0, [qw(A T 12)],   'Chr2:12' ],
    [ $vcf_string, q{1}, 'Sample1', undef, 0, [qw(A T 12)], 'Chr1:12 Sample1' ],
    [ $vcf_string, q{1}, 'Sample1', undef, 1, [qw(G T 14)], 'Chr1:14 Sample1' ],
    [ $vcf_string, q{1}, 'Sample1', undef, 2, [qw(A T 16)], 'Chr1:16 Sample1' ],
    [
        $vcf_string, q{1}, 'Sample2', undef,
        0, [qw(G C,T 14)], 'Chr1:14 Sample2'
    ],
    [ $vcf_string, q{1}, 'Sample2', undef, 1, [qw(A C 16)], 'Chr1:16 Sample2' ],
    [ $vcf_string, q{1}, 'Sample2', undef, 2, [qw(A T 20)], 'Chr1:20 Sample2' ],
);
use warnings;

my @get_variants_by_chr_throws_data = (
    [ undef, q{1}, undef, undef, 'VCF filehandle argument missing' ],
    [
        $vcf_string, undef, undef, undef,
        'Required chromosome argument missing'
    ],
    [ $vcf_string_noheader, q{1}, 'Sample1', undef, 'VCF has no header' ],
    [ $vcf_string,          q{1}, 'SampleX', undef, 'Unknown sample' ],
    [
        $vcf_string,       q{1},
        'SampleX,SampleY', undef,
        'Unknown sample',  'Unknown samples'
    ],
    [
        $vcf_string,       q{1},
        'Sample1,SampleX', undef,
        'Unknown sample',  'Unknown and known sample'
    ],
    [
        $vcf_string_diff_ref, q{1},
        undef,                undef,
        'Different REF at same position'
    ],
);

my @edit_seq_like_data = (
    [ 'AGCT', q{G}, q{C}, 2, $METHOD{ALT}, 'ACCT', 'Simple ALT' ],
    [
        'AGCT',               q{G},      q{C}, 2,
        $METHOD{NOT_REF_ALT}, 'A[AT]CT', 'Simple NOT_REF_ALT'
    ],
    [ 'AgcT', q{G}, q{C}, 2, $METHOD{ALT}, 'AccT', 'Preserve case simple ALT' ],
    [
        'AgcT', q{G}, q{C}, 2, $METHOD{NOT_REF_ALT}, 'A[at]cT',
        'Preserve case simple NOT_REF_ALT'
    ],
    [ 'AGCT', q{G}, [qw(C T)], 2, $METHOD{ALT}, 'A[CT]CT', 'Multiple ALT' ],
    [
        'AGCT', q{G}, [qw(C T)], 2,
        $METHOD{NOT_REF_ALT}, 'AACT', 'Multiple NOT_REF_ALT'
    ],
    [
        'AgCT', q{G}, [qw(C T)], 2, $METHOD{ALT}, 'A[ct]CT',
        'Preserve case multiple ALT'
    ],
    [
        'AgCT', q{G}, [qw(C T)], 2, $METHOD{NOT_REF_ALT}, 'AaCT',
        'Preserve case multiple NOT_REF_ALT'
    ],
    [ 'AGCT', q{G}, 'A,T', 2, $METHOD{ALT}, 'A[AT]CT', 'Comma ALT' ],
    [
        'AGCT',               q{G},   'A,T', 2,
        $METHOD{NOT_REF_ALT}, 'ACCT', 'Comma NOT_REF_ALT'
    ],
    [
        'AgCT', q{G}, 'A,T', 2, $METHOD{ALT}, 'A[at]CT',
        'Preserve case comma ALT'
    ],
    [
        'AgCT', q{G}, 'A,T', 2, $METHOD{NOT_REF_ALT}, 'AcCT',
        'Preserve case comma NOT_REF_ALT'
    ],
    [ 'ANCT', q{N}, q{C}, 2, $METHOD{ALT}, 'ANCT', 'N unchanged ALT' ],
    [
        'ANCT',               q{N},   q{C}, 2,
        $METHOD{NOT_REF_ALT}, 'ANCT', 'N unchanged NOT_REF_ALT'
    ],
    [
        'AncT', q{N}, q{C}, 2, $METHOD{ALT}, 'AncT',
        'Preserve case N unchanged ALT'
    ],
    [
        'AncT', q{N}, q{C}, 2, $METHOD{NOT_REF_ALT}, 'AncT',
        'Preserve case N unchanged NOT_REF_ALT'
    ],
    [ 'AGCT', q{G}, '', 2, $METHOD{ALT}, 'AGCT', 'No ALT' ],
    [
        'AGCT', q{G}, 'A,C,T', 2,
        $METHOD{NOT_REF_ALT}, 'AGCT', 'No NOT_REF_ALT'
    ],
);

my @edit_seq_throws_data = (
    [ undef,  q{G},  q{C},  2,     $METHOD{ALT}, 'Sequence argument missing' ],
    [ 'AGCT', undef, q{C},  2,     $METHOD{ALT}, 'REF argument missing' ],
    [ 'AGCT', q{G},  undef, 2,     $METHOD{ALT}, 'ALT argument missing' ],
    [ 'AGCT', q{G},  q{C},  undef, $METHOD{ALT}, 'Position argument missing' ],
    [ 'AGCT', q{G},  q{C},  2,     undef,        'Method argument missing' ],
    [
        'AGCT', q{G}, q{C}, 5, $METHOD{ALT},
        'Invalid position',
        'Position too large'
    ],
    [
        'AGCT', q{G}, q{C}, 0, $METHOD{ALT},
        'Invalid position',
        'Position too small'
    ],
    [ 'AGCT', q{G}, q{C}, q{x}, $METHOD{ALT}, 'Invalid position' ],
    [
        'AGCT', 'GC', q{C}, 2, $METHOD{ALT},
        'Multiple base REF not yet supported',
        'Multiple base REF'
    ],
    [
        'ARCT', q{R}, q{G}, 2, $METHOD{ALT}, 'REF can only be A, G, C, T or N',
        'Invalid REF'
    ],
    [
        'AGCT', q{C}, q{G}, 2, $METHOD{ALT},
        'REF does not match genome',
        'REF not in genome'
    ],
    [ 'AGCT', q{G}, q{C}, 2, 0, 'Unknown method argument' ],
);

plan tests => ( 2 * scalar @edit_genome_like_data ) +
  scalar @edit_genome_throws_data +
  ( 3 * scalar @get_chr_is_data ) +
  scalar @get_chr_throws_data +
  scalar @get_variants_by_chr_is_data +
  scalar @get_variants_by_chr_is_deeply_data +
  scalar @get_variants_by_chr_throws_data +
  scalar @edit_seq_like_data +
  scalar @edit_seq_throws_data;

foreach my $datum (@edit_genome_like_data) {
    my ( $genome, $vcf, $method, $regex, $name ) = @{$datum};
    my $output;
    open my $fh, q{>}, \$output;
    edit_genome(
        { genome => $genome, vcf => $vcf, method => $method, fh => $fh } );
    like( $output, qr/\A$regex\z/ms, $name . ' filehandle' );
    close $fh;
}

foreach my $datum (@edit_genome_like_data) {
    my ( $genome, $vcf, $method, $regex, $name ) = @{$datum};
    my $output = capture_stdout {
        edit_genome( { genome => $genome, vcf => $vcf, method => $method } );
    };
    like( $output, qr/\A$regex\z/ms, $name . ' STDOUT' );
}

foreach my $datum (@edit_genome_throws_data) {
    my ( $genome, $vcf, $method, $regex, $name ) = @{$datum};
    if ( !defined $name ) {
        $name = $regex;
    }
    throws_ok {
        edit_genome( { genome => $genome, vcf => $vcf, method => $method } )
    }
    qr/$regex/ms, $name;
}

my $prev_genome = q{};
my $genome_fh;
foreach my $datum (@get_chr_is_data) {
    my ( $genome, $short_name, $long_name, $sequence, $name ) = @{$datum};
    if ( $prev_genome ne $genome ) {
        open $genome_fh, q{<}, \$genome;
    }
    $prev_genome = $genome;
    my ( $got_short_name, $got_long_name, $got_sequence ) = get_chr($genome_fh);
    is( $got_short_name, $short_name, $name . ' Short Name' );
    is( $got_long_name,  $long_name,  $name . ' Long Name' );
    is( $got_sequence,   $sequence,   $name . ' Sequence' );
}
close $genome_fh;

foreach my $datum (@get_chr_throws_data) {
    my ( $genome, $regex, $name ) = @{$datum};
    if ( !defined $name ) {
        $name = $regex;
    }
    my $fh;
    open $fh, q{<}, \$genome if $genome;
    throws_ok { get_chr($fh) } qr/$regex/ms, $name;
    close $fh if $genome;
}

foreach my $datum (@get_variants_by_chr_is_data) {
    my ( $vcf, $chr, $samples, $only_pass, $count, $name ) = @{$datum};
    $vcf =~ s/[ ]+/\t/xmsg;    # Ensure VCF is tab separated
    open my $fh, q{<}, \$vcf;
    my $variants = get_variants_by_chr( $fh, $chr, $samples, $only_pass );
    is( scalar @{$variants}, $count, $name );
    close $fh;
}

foreach my $datum (@get_variants_by_chr_is_deeply_data) {
    my ( $vcf, $chr, $samples, $only_pass, $idx, $variant, $name ) = @{$datum};
    $vcf =~ s/[ ]+/\t/xmsg;    # Ensure VCF is tab separated
    open my $fh, q{<}, \$vcf;
    my $variants = get_variants_by_chr( $fh, $chr, $samples, $only_pass );
    is_deeply( $variants->[$idx], $variant, $name );
    close $fh;
}

foreach my $datum (@get_variants_by_chr_throws_data) {
    my ( $vcf, $chr, $samples, $only_pass, $regex, $name ) = @{$datum};
    if ( !defined $name ) {
        $name = $regex;
    }
    my $fh;
    if ($vcf) {
        $vcf =~ s/[ ]+/\t/xmsg;    # Ensure VCF is tab separated
        open $fh, q{<}, \$vcf;
    }
    throws_ok { get_variants_by_chr( $fh, $chr, $samples, $only_pass ) }
    qr/$regex/ms, $name;
    close $fh if $vcf;
}

foreach my $datum (@edit_seq_like_data) {
    my ( $seq, $ref, $alt, $pos, $method, $regex, $name ) = @{$datum};
    my $seq_ref = \$seq;
    edit_seq( $seq_ref, $ref, $alt, $pos, $method );
    like( $seq, qr/\A $regex \z/xms, $name );
}

foreach my $datum (@edit_seq_throws_data) {
    my ( $seq, $ref, $alt, $pos, $method, $regex, $name ) = @{$datum};
    my $seq_ref;
    if ( defined $seq ) {
        $seq_ref = \$seq;
    }
    if ( !defined $name ) {
        $name = $regex;
    }
    throws_ok { edit_seq( $seq_ref, $ref, $alt, $pos, $method ) } qr/$regex/ms,
      $name;
}
