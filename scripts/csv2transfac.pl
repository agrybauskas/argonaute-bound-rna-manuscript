#!/usr/bin/perl

use strict;
use warnings;

use autodie;
use Getopt::Long;

Getopt::Long::Configure qw( gnu_getopt );

#*csv2transfac.pl
#*    Tries to convert nucleotide frequency table from csv to transfac.
#*
#*Usage:
#*    csv2transfac.pl input.csv > output.transfac
#*

# --------------------------------- Options ----------------------------------- #

my $no_header = 0;

GetOptions(
    'skip-header|s' => sub { $no_header = 1 },
    'help|h' => sub {
                        open my $fh, '<', $0;
                        while( <$fh> ) {
                            my ( $help_message ) = $_ =~ m/^#(\*.*)/;
                            if( $help_message ) {$help_message =~ s/^\*$/ /smgx;}
                            if( $help_message ) { $help_message =~ s/^\*//smgx; }
                            if( $help_message ) { print $help_message, "\n"; }
                        }
                        close $fh;
                        exit;
                    },
) or die "mistake in command line argument.\n";

# ----------------------------------- Main ------------------------------------ #

local @ARGV = ( q{-} ) unless @ARGV;

my %frequency = ();
my $row_count = 1;
while( <> ) {
    if( $row_count == 1 && $no_header ) {
        $row_count++;
        next;
    }

    my ( $nucl, $pos, $freq ) = split /[,\t]/, $_;

    $nucl =~ s/\s//g;
    $pos =~ s/\s//g;
    $freq =~ s/\s//g;

    if( ! defined $nucl || ! defined $pos || ! defined $freq ) {
        die "First three columns are mandatory in csv file.";
    }

    $frequency{$pos}{$nucl} = sprintf '%.0f', $freq * 100;

    $row_count++;
}

my ( $any_pos ) = sort keys %frequency;

print "P0\t" . join( "\t", sort keys %{ $frequency{$any_pos} } ) . "\n";

for my $pos ( sort { $a <=> $b } keys %frequency ) {
    print $pos . "\t" . join( "\t",
                              map { $frequency{$pos}{$_} }
                              sort keys %{ $frequency{$pos} } ) . "\n";
}
