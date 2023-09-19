#!/usr/bin/perl
# biogeochemTo.pl 1322 2020-03-17 21:06:11Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
#
#  Print banner.
#
print("\nWelcome to biogeochemTo ...\n\n");
print("This perl program transforms BioGeoChem output into\n");
print("formatted input for equilibrium, conservation, and kinetic\n");
print("equation chemsitry for STOMP.\n");
#
#  Check command line for options.
#
while ( $ARGV[0] =~ /^-/i ) {
  if( $ARGV[0] =~ /help\b/i ) {
    print "Command Line Entry:\n";
    print "    [Options]\n";
#    print "    STOMP Input File Name\n\n";
    print "    BioGeoChem Output File Name\n";
    print "Options:\n\n";
    print "-help      this help\n";
    print "Examples:\n\n";
    print "biogeochemTo.pl biogeochem.out > stomp.input \n";
    die   "  -- input file name will be prompted\n\n";
#
#  Unrecognized option
#
  } else {
    die "Error: Unrecognized Option: $ARGV[0]\n";
  }
}
if( $ARGV[0] ) {
  $in_file = $ARGV[0];
} else {
  do {
    print("BioGeoChem Output File Name?\n");
    $in_file = <STDIN>;
    chomp( $in_file );
    if( -r $in_file ) {
      $stops = 1;
    } elsif( -B $in_file )  {
      $stops = 0;
      print("Error: BioGeoChem Output File is Binary: $in_file[$nf].\n");
    } else {
      $stops = 0;
       print("Error: BioGeoChem Output File is Unreadable: $in_file[$nf].\n");
    }
    if( $stops == 0 ) { print("Try again!\n\n"); }
  } until $stops != 0;
}
#
#  Open BioGeoChem output file
#

open( OUTPUT,$in_file ) || die "Error: Unable to Open BioGeoChem Output File: $in_file.\n";
@output_array = <OUTPUT>;
#
#  Initialize flags and variables
#
$rd_spnm = 0;
$rd_stoi = 0;
$rd_cmsp = 0;
$rd_comp = 0;
$rd_eqsp = 0;
$rd_equi = 0;
$rd_knvr = 0;
$rd_kine = 0;
$neqc = 0;
$neqe = 0;
$neqk = 0;
#
#  Loop over lines in BioGeoChem output files
#
foreach $output_line (@output_array) {
#
#  Remove return from line
#
  chomp( $output_line );
#
#  Remove leading blank spaces from line
#
  $output_line =~ s/^\s+//;
#
#  Read specie names
#
  if( $rd_spnm >= 1 ) {
    @fields = split(/\s+/,$output_line);
    if( $fields[0] ) {
      push(@spnm,$fields[1]);
      push(@spch,$fields[3]);
    } else {
      $rd_spnm = -1;
    }
  }
#
#  Flag species name read
#
  if( $output_line =~ /SP_name/i) {
    $rd_spnm++;
  }
#
#  Read reaction stoichiometric coefficients
#
  if( $rd_stoi >= 5 ) {
    @fields = split(/\s+/,$output_line);
    if( $fields[0] ) {
      $nsp = ($#fields-1)/2;
      push(@stoi,$fields[0]);
      push(@stoi,$nsp);
      for( $i = 2; $i <= $#fields; $i++ ) {
        push(@stoi,$fields[$i]);
      }
    } else {
      $rd_stoi = -1;
    }
  }
#
#  Flag stoichiometric read and skip four lines
#
  if( $rd_stoi >= 1 && $rd_stoi <= 4 ) {
    $rd_stoi++;
  }
  if( $output_line =~ /\*\*\*\*INPUT STOIC/i) {
    $rd_stoi++;
  }
#
#  Read the number of component species
#
  if( $output_line =~ /NO. OF COMPONENT SPECIES:/i) {
    @fields = split(/\s+/,$output_line);
    $nspc = $fields[5];
  }
#
#  Read the number of equilibrium reactions
#
  if( $output_line =~ /NO. OF EQUILIBRIUM SPECIES:/i) {
    @fields = split(/\s+/,$output_line);
    $nspe = $fields[5];
  }
#
#  Read the number of kinetic reactions
#
  if( $output_line =~ /NO. OF KINETIC SPECIES:/i) {
    @fields = split(/\s+/,$output_line);
    $nspk = $fields[5];
  }
#
#  Read component specie constituents
#
  if( $rd_comp >= 4 ) {
    @fields = split(/\s+/,$output_line);
    push(@comp,$fields[1]);    
    push(@comp,$fields[2]);    
    $ncomp--;
    if( $ncomp <= 0 ) {
      $rd_comp = 0;
    }
  }
#
#  Flag component read and skip three lines
#
  if( $rd_comp >= 1 && $rd_comp <= 3 ) {
    $rd_comp++;
  }
  if( $rd_cmsp > 0 && $output_line =~ /TOTAL SPECIES INVOLVED:/i) {
    $rd_comp++;
    @fields = split(/\s+/,$output_line);
    $nf = $#fields;
    push(@comp,$fields[4]);
    push(@comp,$fields[$nf]);
    $ncomp = $fields[$nf];
    $neqc++;
  }
#
#  Flag component species read
#
  if( $output_line =~ /\*\*\*\*COMPONENT SPECIES/i) {
    $rd_cmsp++;
  }
#
#  Read equilibrium specie constituents
#
  if( $rd_equi >= 1 && $rd_equi <= $nesp && $eskip == 0 ) {
    @fields = split(/\s+/,$output_line);
    push(@equi,$fields[1]);    
    push(@equi,$fields[2]);
    $rd_equi++
  }
  if( $rd_eqsp > 0 && $output_line =~ /log OF EQUILIBRIUM CONSTANT:/i) {
    @fields = split(/\s+/,$output_line);
    push(@equi,$fields[4]);
  }
  if( $rd_eqsp > 0 && $output_line =~ /TOT NO. OF SPECIES:/i) {
    $rd_equi = 1;
    @fields = split(/\s+/,$output_line);
    $nf = $#fields;
    push(@equi,$fields[$nf]+1);
    push(@equi,$fields[$nf-5]);
    $nesp = $fields[$nf];
    $neqe++;
    $eskip = 4;
  }
#
#  Flag equilibrium species read
#
  if( $output_line =~ /\*\*\*EQUILIBRIUM SPECIES/i) {
    $rd_eqsp++;
  }
#
#  Read kinetic species
#
  if( $rd_kine >= 1 && $rd_kine <= $nknv && $kskip == 0 ) {
    @fields = split(/\s+/,$output_line);
    push(@kine,$fields[1]);    
    push(@kine,$fields[2]);
    if( $rd_kine == $nknv ) {
      $kskip = 4;
    }
    $rd_kine++;
  }
#
#  Read kinetic reactions
#
  if( $rd_kine >= 1+$nknv && $rd_kine <= $nknv+$nknr && $kskip == 0  ) {
    @fields = split(/\s+/,$output_line);
    push(@kine,$fields[0]);
    $ifind = 0;
    for( $i = 0; $i < $#krgn; $i++ ) {
      if( $krgn[$i] == $fields[0] ) {
        $ifind = 1;
      }
    }
    if( $ifind == 0 ) {
      push(@krgn,$fields[0]);
    }
    push(@kine,$fields[1]);    
    $rd_kine++;
  }
#
#  Read number of kinetic species and number of reactions
#
  if( $rd_knvr > 0 && $output_line =~ /NO. OF VARIABLES IN dE\/dT:/i) {
    $rd_kine = 1;
    @fields = split(/\s+/,$output_line);
    $nf = $#fields;
    push(@kine,$fields[5]);
    push(@kine,$fields[$nf]);
    $nknv = $fields[5];
    $nknr = $fields[$nf];
    $neqk++;
    $kskip = 4;
  }
#
#  Flag kinetic variables and reactions read
#
  if( $output_line =~ /\*\*\*KINETIC VARIABLES AND REACTIONS INVOLVED/i) {
    $rd_knvr++;
  }
#
#  Skip lines
#
  if( $eskip > 0 ) {
    $eskip--;
  }
  if( $kskip > 0 ) {
    $kskip--;
  }
}
#
#  Print Aqueous Species Card
#
if( $#spnm ) {
  print "\n~Aqueous Species Card\n";
  $nsp = $#spnm+1;
  print "$nsp,1.e-9,cm^2/s,Constant Activity,1.0,\n";
  for( $i = 0; $i <= $#spnm; $i++ ) {
    print "$spnm[$i],$spch[$i],0.0,A,0.0,kg/kmol,\n";
  }
}
#
#  Print Convervation Equations Card
#
if( $neqc ) {
  print "\n~Conservation Equations Card\n";
  print "$neqc,\n";
  $i = 0;  
  while( $i <= $#comp ) {
    $icomp = $comp[$i];
    $i++;
    print "Total_$spnm[$comp[$i-1]-1],";
    print "$comp[$i],";
    $nsp = $comp[$i]-1;
    $i++;
#
#  Print Equation Component First
#
    $ix = $i;
    for( $j = 0; $j <= $nsp; $j++ ) {
      if( $icomp == $comp[$i] ) {
        print "$spnm[$comp[$i]-1],";
        $i++;
        print "$comp[$i],";
        $i++;
      } else {
        $i++;
        $i++;
      }
    }
    $i = $ix;
    for( $j = 0; $j <= $nsp; $j++ ) {
      if( $icomp != $comp[$i] ) {
        print "$spnm[$comp[$i]-1],";
        $i++;
        print "$comp[$i],";
        $i++;
      } else {
        $i++;
        $i++;
      }
    }
    print "\n";
  }
}
#
#  Print Equilibrium Cards
#
if( $neqe ) {
#
#  Print Equilibrium Reactions Card
#
  print "\n~Equilibrium Reactions Card\n";
  print "$neqe,\n";
  $i = 0;
  $er = 0;
  while( $i <= $#equi ) {
    $er++;
    print "EqRc-$er,";
    $nsp = $equi[$i]-1;
    for( $j = 0; $j <= $nsp; $j++ ) {
      $i++;
      $i++;
    }
    print "0.0,$equi[$i],0.0,0.0,0.0,1/mol,\n";
    $i++;
  }
#
#  Print Equilibrium Equations Card
#
  print "\n~Equilibrium Equations Card\n";
  print "$neqe,\n";
  $i = 0;
  $er = 0;
  while( $i <= $#equi ) {
    $er++;
    print "$equi[$i],";
    $nsp = $equi[$i]-1;
    $i++;
    print "$spnm[$equi[$i]-1],";
    $i++;
    for( $j = 1; $j <= $nsp; $j++ ) {
      print "$spnm[$equi[$i]-1],";
      $i++;
      print "$equi[$i],";
      $i++;
    }
    print "EqRc-$er,1.0,\n";
    $i++;
  }
}
#
#  Print Kinetic Cards
#
if( $neqk ) {
#
#  Print Kinetic Reactions Card
#
  print "\n~Kinetic Reactions Card\n";
  $nrck = $#krgn+1;
  print "$nrck,\n";
  $i = 0;
  while( $i <= $#krgn ) {
    print "KnRc-$krgn[$i],--- Add Kinetic Reaction ---\n";
    $i++;
  }
#
#  Print Kinetic Equations Card
#
  print "\n~Kinetic Equations Card\n";
  print "$neqk,\n";
  $i = 0;
  while( $i <= $#kine ) {
    print "Kinetic_$spnm[$kine[$i+2]-1],";
    $nsp = $kine[$i];
    print "$nsp,";
    $i++;
    $nrc = $kine[$i];
    $i++;
    for( $j = 1; $j <= $nsp; $j++ ) {
      print "$spnm[$kine[$i]-1],";
      $i++;
      print "$kine[$i],";
      $i++;
    }
    print "\n";
    print "$nrc,";
    for( $j = 1; $j <= $nrc; $j++ ) {
      print "KnRc-$kine[$i],";
      $i++;
      print "$kine[$i],";
      $i++;
    }
    print "\n";
  }
}
