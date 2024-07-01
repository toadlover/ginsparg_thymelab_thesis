#!/usr/bin/perl -w
#=======

#by Sini 200705  
#Sinisa Bjelic, Ph.D.

use Cwd;
use Getopt::Long;
# Find program options

#dafault values
$protchain="XX";
$ligchain="XX";

#command line values 
GetOptions ('pdbfile=s'=>\$pdb,
            'protchain:s'=>\$protchain,
            'ligchain:s'=>\$ligchain,
             );

open (PDB, "< $pdb");
$ligname="NOLOG";
foreach (<PDB>) {
  $line=$_;
  if ($line =~ /^HETNAM\s+(\w+)\s+/) {
   #HETNAM     UP6 6-AZA URIDINE 5'-MONOPHOSPHATE
    $temp=$1;
    #print "$temp\n"; 
    if ( ($temp ne "CL") && ($temp ne "DOD") ){ $ligname=$1; }
  }
  
  if ($line =~ /^ATOM\s+\d+\s+\w+\s+\w+\s+$protchain/) {
   #ATOM    651   CB    ALA A  89
   print "$line";
  }

  if ($line =~ /^(HETATM\s+\d+\s+\w+\s+\w+\s+)$ligchain(\s+\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)/) {
#                 HETATM    1     C1    DM1          X            1       26.211       30.895       26.958        1.00        20.00      C

                  #ATOM    651    CB                      ALA         A          89          21.991       22.367       -3.427        1.00        37.75           O
  #print "$line";
    print "$1X$2\n";
  }

}
