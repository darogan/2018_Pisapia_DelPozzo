#!/usr/bin/perl
#------------------------------------------------------------------------------
# Perl Data Analysis:
# Pisapia et al (2018)
# Tristetraprolin regulates the turnover of autoimmune-associated HLA-DQ mRNAs 
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/darogan/2018_Pisapia_DelPozzo
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics, 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

use strict;


my %Sequences;

# WT
$Sequences{"DQA101"} = "UGAAUCCCAUCCUGGAAGGGAAGUGCAUCGCCAUCUACAGGAGCAGAAGAGUGGACUUGCUACAUGACCUAGCACUAUUCUCUGGCCCGAUUUAUCAUAUCCCUUUUCUCCUCCAAAUAUUUCUCCUCUCACCUUUUCUGUGGGACUUAAGCUGCUAUAUCCCCUCAGAGCUCACAAAUGUCUUU";
$Sequences{"DQB105"} = "UGACUCCUGAGACUGUUUUAACUAAGACUGGUUAUCACUCUUCUGUGAUGCCUGCUUGUCCCUGCCCAGAAUUCCCAGCUGCCUGUGUCAGCUUGUCCCCCUGAGAUCAAAGUCCUACAGUGGCUGUCACGCAACCACCAGGUCAUCUCCUUUCAUCCCCACCCCAAGGCGCUGGCUGUGACUCUGCUUCCUGCACUGACCCAGAGCC";

# CD
$Sequences{"DQA105"} = "UGAAUCCCAUCCUGGAAUGGAAGUGCAUCGCCAUCUACAGGAGCAGAAGAGUGGACUUGCUACAUGACCUAGCAUUAUUUUCUGGCCCCAUUUAUCAUAUCCCUUUUCUCCUCCAAAUGUUUCUCCUCUCACCUCUUCUGUGGGACUUAAAUUGCUAUAUCUGCUCAGAGCUCACAAAUGCCUUU";
$Sequences{"DQB102"} = "UGACUCCUGAGACUAUUUUAACUGGGAUUGGUUAUCACUUUUCUGUAACGCCUGCUUGUCCCUGCCCAGAAUUCCCAGCUGUCUGUGUCAGCCUGUCCCCCUGAGAUCAGAGUCCUACAGUGGCUGUCACGCAGCCACCAGGUCAUCUCCUUUCAUCCCCACCUUGAGGCGGAUGGCUGUGACCCUACUUCCUGCACUGACCCACAGCC";

my ($i, $genome,
    $dihead, $ditail, $trihead, $tritail, $quadhead, $quadtail, 
    $pentahead, $pentatail, $hexahead, $hexatail,
    $mononucFile, $dinucFile, $trinucFile, $quadnucFile, $pentanucFile, $hexanucFile,
    %mono_nt, %di_nt, %tri_nt, %quad_nt, %penta_nt, %hexa_nt);

$mononucFile  = "mononuc.table.txt";
$dinucFile    = "dinuc.table.txt";
$trinucFile   = "trinuc.table.txt";
$quadnucFile  = "quadnuc.table.txt";
$pentanucFile = "pentanuc.table.txt";
$hexanucFile  = "hexanuc.table.txt";

open(MONO, ">$mononucFile")  || die "Can't open $mononucFile for writing: $!\n";
open(DI,   ">$dinucFile")    || die "Can't open $dinucFile for writing: $!\n";
open(TRI,  ">$trinucFile")   || die "Can't open $trinucFile for writing: $!\n";
open(QUAD, ">$quadnucFile")  || die "Can't open $quadnucFile for writing: $!\n";
open(PENTA,">$pentanucFile") || die "Can't open $pentanucFile for writing: $!\n";
open(HEXA, ">$hexanucFile")  || die "Can't open $hexanucFile for writing: $!\n";

print MONO  "Sequence\tMononucleotide\tObs_freq\tExp_freq\n";
print DI    "Sequence\tDinucleotide\tObs_freq\tExp_freq\n";
print TRI   "Sequence\tTrinucleotide\tObs_freq\tExp_freq\n";
print QUAD  "Sequence\tQuadnucleotide\tObs_freq\tExp_freq\n";
print PENTA "Sequence\tPentanucleotide\tObs_freq\tExp_freq\n";
print HEXA  "Sequence\tHexanucleotide\tObs_freq\tExp_freq\n";

foreach my $seqheader (keys %Sequences)
  {
    print $seqheader, " -- ", $Sequences{$seqheader}, "\n";
	
    my $seq = $Sequences{$seqheader};
    $genome .= uc $seq;
    $dihead = uc substr($seq, 0, 1);
    $di_nt{"$ditail$dihead"}-- if $ditail;
    $ditail = uc substr($seq, -1);

    $trihead = uc substr($seq, 0, 2);
    $tri_nt{"$tritail$trihead"}-- if $tritail;
    $tritail = uc substr($seq, -1);

    $quadhead = uc substr($seq, 0, 3);
    $quad_nt{"$quadtail$quadhead"}-- if $quadtail;
    $quadtail = uc substr($seq, -1);

    $pentahead = uc substr($seq, 0, 4);
    $penta_nt{"$pentatail$pentahead"}-- if $pentatail;
    $pentatail = uc substr($seq, -1);

    $hexahead = uc substr($seq, 0, 5);
    $hexa_nt{"$hexatail$hexahead"}-- if $hexatail;
    $hexatail = uc substr($seq, -1);

    my $len = length $genome;
    for my $i (0..$len-2) 
      {
        my $each_mono_nt = substr($genome, $i, 1);
        my $each_di_nt   = substr($genome, $i, 2);
        $mono_nt{$each_mono_nt}++;
        $di_nt{$each_di_nt}++;
      }

    for my $i (0..$len-3)
      {
        my $each_tri_nt  = substr($genome, $i, 3);
        $tri_nt{$each_tri_nt}++;
      }

    for my $i (0..$len-4)
      {
        my $each_quad_nt  = substr($genome, $i, 4);
        $quad_nt{$each_quad_nt}++;
      }

    for my $i (0..$len-5)
      {
        my $each_penta_nt  = substr($genome, $i, 5);
        $penta_nt{$each_penta_nt}++;
      }

    for my $i (0..$len-6)
      {
        my $each_hexa_nt  = substr($genome, $i, 6);
        $hexa_nt{$each_hexa_nt}++;
      }

    $mono_nt{$ditail}++;

    print "+ Calculating Single nucleotide frequency...\n";
    for my $nt (sort keys %mono_nt) 
      {
        print MONO "$seqheader\t$nt\t", $mono_nt{$nt} / $len, "\t0.25\n";
      }

    print "+ Calculating Dinucleotide frequency...\n";
    for my $nt_pair (sort keys %di_nt) 
      {
        my ($first_nt, $second_nt) = split //, $nt_pair;
        print DI "$seqheader\t$nt_pair\t", $di_nt{$nt_pair} / ($len-1), "\t",
        $mono_nt{$first_nt} * $mono_nt{$second_nt} /$len /$len, "\n";
      }

    print "+ Calculating Trinucleotide frequency...\n";
    for my $nt_tri (sort keys %tri_nt)
      {
        my ($first_nt, $second_nt, $third_nt) = split //, $nt_tri;
        print TRI "$seqheader\t$nt_tri\t", $tri_nt{$nt_tri} / ($len-1), "\t",
        $mono_nt{$first_nt} * $mono_nt{$second_nt} * $mono_nt{$third_nt} /$len /$len /$len, "\n";
      }

    print "+ Calculating Quadnucleotide frequency...\n";
    for my $nt_quad (sort keys %quad_nt)
      {
        my ($first_nt, $second_nt, $third_nt, $forth_nt) = split //, $nt_quad;
        print QUAD "$seqheader\t$nt_quad\t", $quad_nt{$nt_quad} / ($len-1), "\t",
              $mono_nt{$first_nt} * $mono_nt{$second_nt} * $mono_nt{$third_nt} * 
              $mono_nt{$forth_nt} /$len /$len /$len /$len, "\n";
      }

    print "+ Calcualting Pentanucleotide frequency...\n";
    for my $nt_penta (sort keys %penta_nt)
      {
        my ($first_nt, $second_nt, $third_nt, $forth_nt, $fifth_nt) = split //, $nt_penta;
        print PENTA "$seqheader\t$nt_penta\t", $penta_nt{$nt_penta} / ($len-1), "\t",
              $mono_nt{$first_nt} * $mono_nt{$second_nt} * $mono_nt{$third_nt} * 
              $mono_nt{$forth_nt} * $mono_nt{$fifth_nt} /$len /$len /$len /$len /$len, "\n";
      }

    print "+ Calcualting Hexanucleotide frequency...\n";
    for my $nt_hexa (sort keys %hexa_nt)
      {
        my ($first_nt, $second_nt, $third_nt, $forth_nt, $fifth_nt, $sixth_nt) = split //, $nt_hexa;
        print HEXA "$seqheader\t$nt_hexa\t", $hexa_nt{$nt_hexa} / ($len-1), "\t",
              $mono_nt{$first_nt} * $mono_nt{$second_nt} * $mono_nt{$third_nt} *
              $mono_nt{$forth_nt} * $mono_nt{$fifth_nt} * $mono_nt{$sixth_nt} /$len /$len /$len /$len /$len /$len, "\n";
      } 
  }

close DI;
close TRI;
close QUAD;
close PENTA;
close HEXA;

exit;



print "+"x100, "\n";
print "+ Density Matrix\n";
print "+"x100, "\n";

print "+ Creating Density Map\n";
my %SeqDensity;
foreach my $seqheader (keys %Sequences)
  {
    my(@Seq_Density_Di, @Seq_Density_Tri, @Seq_Density_Quad, @Seq_Density_Penta, @Seq_Density_Hexa);
    push @Seq_Density_Di,    (0) x (length($Sequences{$seqheader}));
    push @Seq_Density_Tri,   (0) x (length($Sequences{$seqheader}));
    push @Seq_Density_Quad,  (0) x (length($Sequences{$seqheader}));
    push @Seq_Density_Penta, (0) x (length($Sequences{$seqheader}));
    push @Seq_Density_Hexa,  (0) x (length($Sequences{$seqheader}));

    print $seqheader, "...", "\n";

    my $dipatternstring    = "UA|AU|UU";
    my @dipattern          = split(/\|/, $dipatternstring);

    my $tripatternstring   = "[UA][UA][UA]"; 
    my @tripattern         = split(/\|/, $tripatternstring);

    my $quadpatternstring  = "[UA][UA][UA][UA]";
    my @quadpattern        = split(/\|/, $quadpatternstring);

    my $pentapatternstring = "[UA][UA][UA][UA][UA]";
    my @pentapattern       = split(/\|/, $pentapatternstring);

    my $hexapatternstring = "[UA][UA][UA][UA][UA][UA]";
    my @hexapattern       = split(/\|/, $hexapatternstring);


    #print "Dinucleotide Densities\n";
    foreach my $dipat (@dipattern)
       {
         while($Sequences{$seqheader} =~ m/($dipat)/g) 
           {
             $Seq_Density_Di[pos($Sequences{$seqheader})] += 1;
             $Seq_Density_Di[pos($Sequences{$seqheader})+1] += 1;
           }
       }

    #print "Trinucleotide Densities\n";
    foreach my $tripat (@tripattern)
       { 
         while($Sequences{$seqheader} =~ m/($tripat)/g) 
           {
             my $As = () = $1 =~ /A/gi;
             if($As < length($tripat)/2)
             {
             $Seq_Density_Tri[pos($Sequences{$seqheader})]   += 1;
             $Seq_Density_Tri[pos($Sequences{$seqheader})+1] += 1;
             $Seq_Density_Tri[pos($Sequences{$seqheader})+2] += 1;
             }
           }
       }

    #print "Quadnucleotide Densities\n";
    foreach my $quadpat (@quadpattern)
       { 
         while($Sequences{$seqheader} =~ m/($quadpat)/g) 
           {
             my $As = () = $1 =~ /A/gi;
             if($As < length($quadpat)/2)
             {
             $Seq_Density_Quad[pos($Sequences{$seqheader})]   += 1;
             $Seq_Density_Quad[pos($Sequences{$seqheader})+1] += 1;
             $Seq_Density_Quad[pos($Sequences{$seqheader})+2] += 1;
             $Seq_Density_Quad[pos($Sequences{$seqheader})+3] += 1;
             }
           }
       }

    #print "Pentanucleotide Densities\n";
    foreach my $pentapat (@pentapattern)
       {
         while($Sequences{$seqheader} =~ m/($pentapat)/g)
           {
             my $As = () = $1 =~ /A/gi;
             if($As < length($pentapat)/2)
             {
             $Seq_Density_Penta[pos($Sequences{$seqheader})]   += 1;
             $Seq_Density_Penta[pos($Sequences{$seqheader})+1] += 1;
             $Seq_Density_Penta[pos($Sequences{$seqheader})+2] += 1;
             $Seq_Density_Penta[pos($Sequences{$seqheader})+3] += 1;
             $Seq_Density_Penta[pos($Sequences{$seqheader})+4] += 1;
             }
           }
       }

    #print "Pentanucleotide Densities\n";
    foreach my $hexapat (@hexapattern)
       {
         while($Sequences{$seqheader} =~ m/($hexapat)/g)
           {
             my $As = () = $1 =~ /A/gi;
             if($As < length($hexapat)/2)
             {
             $Seq_Density_Hexa[pos($Sequences{$seqheader})]   += 1;
             $Seq_Density_Hexa[pos($Sequences{$seqheader})+1] += 1;
             $Seq_Density_Hexa[pos($Sequences{$seqheader})+2] += 1;
             $Seq_Density_Hexa[pos($Sequences{$seqheader})+3] += 1;
             $Seq_Density_Hexa[pos($Sequences{$seqheader})+4] += 1;
             $Seq_Density_Hexa[pos($Sequences{$seqheader})+5] += 1;
             }
           }
       }

$SeqDensity{$seqheader}{"di"}    = [@Seq_Density_Di];
$SeqDensity{$seqheader}{"tri"}   = [@Seq_Density_Tri];
$SeqDensity{$seqheader}{"quad"}  = [@Seq_Density_Quad];
$SeqDensity{$seqheader}{"penta"} = [@Seq_Density_Penta];
$SeqDensity{$seqheader}{"hexa"}  = [@Seq_Density_Hexa];
}

my @motif_ranges = ("di", "tri", "quad", "penta", "hexa");

print "+ Writing DQA101_DQA105_densities.txt\n";
open(PAIR1,">DQA101_DQA105_densities.txt") || die "Can't open DQA101_DQA105_densities.txt for writing\n";

my @Headers      = ("DQA101", "DQA105");
     
print PAIR1 "position,";
for(my $k=0; $k<=$#Headers; $k++)
   { 
     for(my $i=0; $i<=$#motif_ranges; $i++)
        {
          print PAIR1 "${Headers[$k]}_${motif_ranges[$i]}";
          print PAIR1 ",", if( $k + $i < $#Headers + $#motif_ranges);
        }
   }
print PAIR1 "\n";     

for(my $j=0; $j<=$#{ $SeqDensity{$Headers[0]}{$motif_ranges[$i]} }; $j++)
   {  
     for(my $k=0; $k<=$#Headers; $k++)
        { 
          for(my $i=0; $i<=$#motif_ranges; $i++)
             {
                print PAIR1 "$j,", if(($k == 0) and ($i ==0));
                print PAIR1 $SeqDensity{$Headers[$k]}{$motif_ranges[$i]}[$j];
                print PAIR1 ",", if( $k + $i < $#Headers + $#motif_ranges); 
             }
        } 
     print PAIR1 "\n";
   }

close PAIR1;


print "+ Writing DQB105_DQB102_densities.txt\n";
open(PAIR2,">DQB105_DQB102_densities.txt") || die "Can't open DQA101_DQA105_densities.txt for writing\n";

my @Headers      = ("DQB105", "DQB102");

print PAIR2 "position,";
for(my $k=0; $k<=$#Headers; $k++)
   {
     for(my $i=0; $i<=$#motif_ranges; $i++)
        {
          print PAIR2 "${Headers[$k]}_${motif_ranges[$i]}";
          print PAIR2 ",", if( $k + $i < $#Headers + $#motif_ranges);
        }
   }
print PAIR2 "\n";

for(my $j=0; $j<=$#{ $SeqDensity{$Headers[0]}{$motif_ranges[$i]} }; $j++)
   {
     for(my $k=0; $k<=$#Headers; $k++)
        {
          for(my $i=0; $i<=$#motif_ranges; $i++)
             {
                print PAIR2 "$j,", if(($k == 0) and ($i ==0));
                print PAIR2 $SeqDensity{$Headers[$k]}{$motif_ranges[$i]}[$j];
                print PAIR2 ",", if( $k + $i < $#Headers + $#motif_ranges);
             }
        }
     print PAIR2 "\n";
   }

close PAIR2;

print "+"x100, "\n";
print "+ End of Script\n";
print "+"x100, "\n";
