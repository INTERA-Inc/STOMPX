#!/usr/bin/perl
# plot_bhTo.pl 1080 2017-03-14 16:22:02Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
print("\nWelcome to plot_bhTo ...\n\n");
print("This perl program transforms STOMP plot_bh file(s) into\n");
print("formatted input files for Surfer, Tecplot, and VTK.\n");
print("VTK format can be read by VisIt and Paraview.\n");
print("Multiple plot files are allowable for Tecplot;\n");
print("whereas, Surfer and VTK only accept a single plot file.\n");
print("Entries not made on the command line will be prompted.\n\n");
$nargv = $#ARGV;
#
#  Check command line for options.
#
$t_opt = 0;
$ts_opt = 0;
$ta_opt = 0;
while ( $ARGV[0] =~ /^-/i ) {
  if( $ARGV[0] =~ /help\b/i ) {
    print "Command Line Entry:\n\n";
    print "    Plotting Package Option [Surfer|Tecplot|Tecplot10|VTK]\n";
    print "    PLotting Package File to Generate\n";
    print "    STOMP Borehole Plot File Name(s)\n\n";
    print "Options:\n\n";
    print "-help      this help\n";
    print "-t         prompt for Tecplot title\n";
    print "-ts        use time stamp Tecplot zone names\n";
    print "Plotting Package Options:\n\n";
    print "Surfer\n";
    print "Tecplot\n";
    print "Tecplot10 (for Tecplot Version 10 and older)\n";
    print "VTK\n";
    print "Examples:\n\n";
    print "plot_bhTo.pl -t -ts Tecplot plots_bh.dat plot_bh.123 plot_bh.456 plot_bh.789\n";
    print "plot_bhTo.pl Surfer plot_bh.dat plot_bh.123\n";
    print "plot_bhTo.pl Tecplot plots_bh.dat plot_bh.*\n";
    die "plot_bhTo.pl -t Tecplot plots_bh.dat `ls -1rt plot_bh.*`\n\n";
#
#  Tecplot title option
#
  } elsif( $ARGV[0] =~ /\-t\b/ ) {
    $t_opt = 1;
    shift(@ARGV);
#
#  Tecplot time adjust option
#
  } elsif( $ARGV[0] =~ /\-t\b/ ) {
    $ta_opt = 1;
    shift(@ARGV);
#
#  Tecplot time stamp option
#
  } elsif( $ARGV[0] =~ /\-ts\b/ ) {
    $ts_opt = 1;
    shift(@ARGV);
#
#  Unrecognized option
#
  } else {
    die "Error: Unrecognized Option: $ARGV[0]\n";
  }
}
#
#  Search for plotting package name as first argument or prompt user.
#
if( $ARGV[0] ) {
  $plot_package = $ARGV[0];
#  chomp( $plot_package );
  if( $plot_package =~ /\n/ ) {
    chop( $plot_package );
    print "$plot_package\n";
  }
  if( $plot_package =~ /^tecplot\b/i ) {
    print("Plotting Package: Tecplot\n");
  } elsif( $plot_package =~ /^tecplot10\b/i ) {
    print("Plotting Package: Tecplot10\n");
  } elsif( $plot_package =~ /^surfer\b/i ) {
    print("Plotting Package: Surfer\n");
  } elsif( $plot_package =~ /^vtk\b/i ) {
    print("Plotting Package: VTK\n");
  } else {
    die "Error: Unrecognized Plotting Package: $plot_package.\n\n";
  }
#
#  Search for plotting package file name (output file) as second argument or prompt user.
#
  if( $ARGV[1] ) {
    $out_file = $ARGV[1];
    chomp( $out_file );
    open( OUT,">$out_file") || die "Error: Unable to Create Plotting Package File: $out_file.\n";
    print( "Created Plotting Package File: $out_file.\n");
#
#  Search for plot file name(s) as third(+) argument or prompt user.
#
    if( $ARGV[2] ) {
      print( "Borehole Plot File(s): \n" );
#
#  Check for multiple plot files with the Surfer and VTK plotting packages
#
      if( $plot_package =~ /^surfer\b/i && $#ARGV > 2 ) {
        die "Error: Multiple Borehole Plot Files Specified for the Surfer Plotting Package.\n";
      }
      if( $plot_package =~ /^vtk\b/i && $#ARGV > 2 ) {
        die "Error: Multiple Borehole Plot Files Specified for the VTK Plotting Package.\n";
      }
      for( $nf = 2; $nf <= $#ARGV; $nf++ ) {
        $plot_file[$nf-2] = $ARGV[$nf];
        $nff = $nf-1;
        print( "STOMP Borehole Plot File #$nff = $plot_file[$nf-2]\n");
      }
#
#  No third(+) argument(s); ask user for plot file name(s).
#
    } else {
      do {
        $stops = 1;
        print("STOMP Borehole Plot File Name(s)?\n");
        $names = <STDIN>;
        chomp( $names );
        @plot_file = split(/\s+/,$names);
        $names = join(",",@plot_file);
        @plot_file = split(/:+/,$names);
        $names = join(",",@plot_file);
        @plot_file = split(/,+/,$names);
#
#  Check for multiple plot files with the Surfer or VTK plotting packages
#
        if( $plot_package =~ /^surfer\b/i && $#plot_file > 0 ) {
          die "Error: Multiple Borehole Plot Files Specified for the Surfer Plotting Package.\n";
        }
        if( $plot_package =~ /^vtk\b/i && $#plot_file > 0 ) {
          die "Error: Multiple Borehole Plot Files Specified for the VTK Plotting Package.\n";
        }
        print( "STOMP Borehole Plot File(s): \n" );
        for( $nf = 0; $nf <= $#plot_file; $nf++ ) {
          $nff = $nf+1;
          if( -r $plot_file[$nf] ) {
            print("STOMP Borehole Plot File #$nff = $plot_file[$nf]\n");
          } elsif( -B $plot_file[$nf] )  {
            $stops = 0;
            print("Error: STOMP Borehole Plot File is Binary: $plot_file[$nf].\n");
          } else {
            $stops = 0;
            print("Error: STOMP Borehole Plot File is Unreadable: $plot_file[$nf].\n");
          }
        }
        if( $stops == 0 ) { print("Try again!\n\n"); }
      } until $stops != 0;
    }
#
#  No second argument; ask user for plotting package file name (output file).
#
  } else {
    $stops = 0;
    do {
      if( $plot_package =~ /^tecplot\b/i ) {
        print("Tecplot File to Generate?\n");
      } elsif( $plot_package =~ /^tecplot10\b/i ) {
        print("Tecplot10 File to Generate?\n");
      } elsif( $plot_package =~ /^surfer\b/i ) {
        print("Surfer File to Generate?\n");
      } elsif( $plot_package =~ /^vtk\b/i ) {
        print("VTK File to Generate\n");
      }
      $out_file = <STDIN>;
      chomp( $out_file );
      if( open( OUT,">$out_file") ) {
        $stops = 1;
        print("Plotting Package File: $out_file\n");
      } else {
        print("Error: Unable to Create Plotting Package File: $out_file. Try again!\n\n");
      }
    } until $stops != 0;
#
#  No second or third(+) argument(s); ask user for plot file name(s).
#
    do {
      $stops = 1;
      print("STOMP Borehole Plot File Name(s)?\n");
      $names = <STDIN>;
      chomp( $names );
      @plot_file = split(/\s+/,$names);
      $names = join(",",@plot_file);
      @plot_file = split(/:+/,$names);
      $names = join(",",@plot_file);
      @plot_file = split(/,+/,$names);
#
#  Check for multiple plot files with the Surfer or VTK plotting package
#
      if( $plot_package =~ /^surfer\b/i && $#plot_file > 0 ) {
        die "Error: Multiple Borehole Plot Files Specified for the Surfer Plotting Package.\n";
      }
      if( $plot_package =~ /^vtk\b/i && $#plot_file > 0 ) {
        die "Error: Multiple Borehole Plot Files Specified for the VTK Plotting Package.\n";
      }
      print( "Borehole Plot File(s): \n" );
      for( $nf = 0; $nf <= $#plot_file; $nf++ ) {
        $nff = $nf+1;
        if( -r $plot_file[$nf] ) {
          print("Borehole Plot File #$nff = $plot_file[$nf]\n");
        } elsif( -B $plot_file[$nf] )  {
          $stops = 0;
          print("Error: Borehole Plot File is Binary: $plot_file[$nf].\n");
        } else {
          $stops = 0;
          print("Error: Borehole Plot File is Unreadable: $plot_file[$nf].\n");
        }
      }
      if( $stops == 0 ) { print("Try again!\n\n"); }
    } until $stops != 0;
  }
#
#  No first argument; ask user for plotting package name.
#
} else {
  $stops = 0;
  do {
    print("Plotting Package [Surfer|Tecplot|Tecplot10|VTK]?\n");
    $plot_package = <STDIN>;
    chomp( $plot_package );
    if( $plot_package =~ /^tecplot\b/i ) {
      print("Plotting Package: Tecplot\n");
      $stops = 1;
    } elsif( $plot_package =~ /^tecplot10\b/i ) {
      print("Plotting Package: Tecplot10\n");
      $stops = 1;
    } elsif( $plot_package =~ /^surfer\b/i ) {
      print("Plotting Package: Surfer\n");
      $stops = 1;
    } elsif( $plot_package =~ /^vtk\b/i ) {
      print("Plotting Package: VTK\n");
      $stops = 1;
    } else {
      print("Error: Unrecognized Plotting Package: $plot_package.  Try again!\n\n");
    }
  } until $stops != 0;
#
#  No first or second argument; ask user for output file name.
#
  $stops = 0;
  do {
      if( $plot_package =~ /^tecplot\b/i ) {
        print("Tecplot File to Generate?\n");
      } elsif( $plot_package =~ /^tecplot10\b/i ) {
        print("Tecplot10 File to Generate?\n");
      } elsif( $plot_package =~ /^surfer\b/i ) {
        print("Surfer File to Generate?\n");
      } elsif( $plot_package =~ /^vtk\b/i ) {
        print("VTK File to Generate\n");
      }
    $out_file = <STDIN>;
    chomp( $out_file );
    if( open( OUT,">$out_file") ) {
      $stops = 1;
      print("PLotting Package File to Generate: $out_file\n");
    } else {
      print("Error: Unable to Create Plotting Package File: $out_file. Try again!\n\n");
    }
  } until $stops != 0;
#
#  No first, second or third(+) argument(s); ask user for plot file name(s).
#
  do {
    $stops = 1;
    print("STOMP Borehole Plot File Name(s)?\n");
    $names = <STDIN>;
    chomp( $names );
    @plot_file = split(/\s+/,$names);
    $names = join(",",@plot_file);
    @plot_file = split(/:+/,$names);
    $names = join(",",@plot_file);
    @plot_file = split(/,+/,$names);
#
#  Check for multiple plot files with the Surfer or VTK plotting packages
#
    if( $plot_package =~ /^surfer\b/i && $#plot_file > 0 ) {
      die "Error: Multiple Borehole Plot Files Specified for the Surfer Plotting Package.\n";
    }
    if( $plot_package =~ /^vtk\b/i && $#plot_file > 0 ) {
      die "Error: Multiple Borehole Plot Files Specified for the VTK Plotting Package.\n";
    }
    print( "Borehole Plot File(s): \n" );
    for( $nf = 0; $nf <= $#plot_file; $nf++ ) {
      $nff = $nf+1;
      if( -r $plot_file[$nf] ) {
        print("STOMP Borehole Plot File #$nff = $plot_file[$nf]\n");
      } elsif( -B $plot_file[$nf] )  {
        $stops = 0;
        print("Error: STOMP Borehole Plot File is Binary: $plot_file[$nf].\n");
      } else {
        $stops = 0;
        print("Error: STOMP Borehole Plot File is Unreadable: $plot_file[$nf].\n");
      }
    }
    if( $stops == 0 ) { print("Try again!\n\n"); }
  } until $stops != 0;
}
#
#  Tecplot Title
#
  if( ($plot_package =~ /^tecplot\b/i || $plot_package =~ /^tecplot10\b/i ||
     $plot_package =~ /^tecplotxy\b/i ) 
    && $t_opt == 1 ) {
    print("Tecplot Title?\n");
    $title = <STDIN>;
    chomp( $title );
  }
#
#  Tecplot Time Stamp
#
  if( ($plot_package =~ /^tecplot\b/i || $plot_package =~ /^tecplot10\b/i )
    && $ts_opt == 1 ) {
    $stops = 0;
    do {
      print("Tecplot Time Stamp Unit?\n");
      print("Options: s[econd] | m[inute] | h[our] | d[ay] | w[eek] | y[ear]\n");
      $tsunit = <STDIN>;
      chomp( $tsunit );
      if( $tsunit =~ /^s/i ) {
        $tsunit = "s";
        $stops = 1;
      } elsif( $tsunit =~ /^m/i ) {
        $tsunit = "min";
        $stops = 1;
      } elsif( $tsunit =~ /^h/i ) {
        $tsunit = "h";
        $stops = 1;
      } elsif( $tsunit =~ /^d/i ) {
        $tsunit = "day";
        $stops = 1;
      } elsif( $tsunit =~ /^w/i ) {
        $tsunit = "wk";
        $stops = 1;
      } elsif( $tsunit =~ /^y/i ) {
        $tsunit = "yr";
        $stops = 1;
      } else {
        print("Error: Unrecognized Time Stamp Unit: $tsunit.  Try again!\n\n");
      }
    } until $stops != 0;
  }
#
#  Tecplot Time Adjustment
#
  if( ($plot_package =~ /^tecplot\b/i || $plot_package =~ /^tecplot10\b/i )
    && $ta_opt == 1 ) {
    print("Time Adjustment?\n");
    $tadj = <STDIN>;
    chomp( $tadj );
  }
#
#  Loop over plot files
#
for( $nf = 0; $nf <= $#plot_file; $nf++ ) {
  open( PLOT,$plot_file[$nf] ) || die "Error: Unable to Open STOMP Borehole Plot File: $plot_file[$nf].\n";
  print("Converting STOMP Borehole Plot File: $plot_file[$nf].\n");
  @plot_array = <PLOT>;
#
#  Initialize flags
#
  $xflag = 0;
  $yflag = 0;
  $zflag = 0;
  $nvflag = 0;
  $fvflag = -1;
  $nfv = 0;
  $nvert = 1;
#
#  Loop over lines in plot file
#
  foreach $plot_line (@plot_array) {
#
#  Remove return from line
#
    chomp( $plot_line );
#
#  Remove leading blank spaces from line
#
    $plot_line =~ s/^\s+//;
#
#  Set the number of Borehole Nodes
#
    if( $plot_line =~ /Number of Borehole Nodes/ ) {
      @fields = split(/\s+/,$plot_line);
      $nbn = $fields[$#fields];
#
#  Set the number of vertices - always 8 since these are hexahedra
#
      $nvert = 8;
#
#  Time stamp
#
    } elsif( $plot_line =~ /Time = / && $ts_opt == 1 ) {
      @fields = split(/\s+|,/,$plot_line);
      for( $i = 3; $i <= $#fields; $i++ ) {
        if( $fields[$i] =~ /$tsunit/ ) {
          push(@tsv,$fields[$i-1]);
        }
      }
#
#  Located a character string starting a line
#
    } elsif( $plot_line =~ /^[A-Z]|^[a-z]/ ) {
#
#  Set flag to read x-direction borehole node positions
#
       if( $plot_line =~ /X-Direction Borehole Node Vertices/ ) {
        @fields = split(/,/,$plot_line);
        $xunit = $fields[1];
        $xunit =~ s/\s+//g;
        $xvert = 1;
        $xflag = 1;
#
#  Set flag to read y-direction borehole node positions
#
      } elsif( $plot_line =~ /Y-Direction Borehole Node Vertices/ ) {
        @fields = split(/,/,$plot_line);
        $yunit = $fields[1];
        $yunit =~ s/\s+//g;
        $yvert = 1;
        $yflag = 1;
#
#  Set flag to read z-direction borehole node positions
#
      } elsif( $plot_line =~ /Z-Direction Borehole Node Vertices/ ) {
        @fields = split(/,/,$plot_line);
        $zunit = $fields[1];
        $zunit =~ s/\s+//g;
        $zvert = 1;
        $zflag = 1;
#
#  Set flag to read node volumes
#
      } elsif( $plot_line =~ /Borehole Node Volume/ ) {
        @fields = split(/,/,$plot_line);
        $nvunit = $fields[1];
        $nvunit =~ s/\s+//g;
        $nvflag = 1;
#
#  Set flag to read field variables
#
      } elsif( $fvflag >= 0 ) {
        @fields = split(/,/,$plot_line);
        push(@fvname,$fields[0]);
        $unit = $fields[1];
        $unit =~ s/\s+//g;
        push(@fvunit,$unit);
        $fvflag = 1;
        $nfv++;
      }
#
#  Read x-direction vertex positions
#
    } elsif( $xflag != 0 ) {
      if( $xvert ) {
          @fields = split(/\s+/,$plot_line);
          push(@xp,@fields);
          if( ($#xp+1) >= $nvert*$nbn ){ $xflag = 0; }
      }   else {
          @fields = split(/\s+/,$plot_line);
          push(@xp,@fields);
          if( ($#xp+1) >= $nbn ){ $xflag = 0; }
      }
#
#  Read y-direction vertex positions
#
    } elsif( $yflag != 0 ) {
      if( $yvert ) {
          @fields = split(/\s+/,$plot_line);
          push(@yp,@fields);
          if( ($#yp+1) >= $nvert*$nbn ){ $yflag = 0; }
      }   else {
          @fields = split(/\s+/,$plot_line);
          push(@yp,@fields);
          if( ($#yp+1) >= $nbn ) { $yflag = 0; }
      }
#
#  Read z-direction vertex positions
#
    } elsif( $zflag != 0 ) {
      if( $zvert ) {
          @fields = split(/\s+/,$plot_line);
          push(@zp,@fields);
          if( ($#zp+1) >= $nvert*$nbn ){ $zflag = 0; }
      }   else {
          @fields = split(/\s+/,$plot_line);
          push(@zp,@fields);
          if( ($#zp+1) >= $nbn ) { $zflag = 0; }
      }
#
#  Read Fracture Borehole Node Volume
#
    } elsif( $nvflag != 0 ) {
      @fields = split(/\s+/,$plot_line);
      push(@nv,@fields);
      if( ($#nv+1) >= $nbn ) {
        $nvflag = 0;
        $fvflag = 0;
      }
#
#  Read field variable
#
    } elsif( $fvflag > 0 ) {
      @fields = split(/\s+/,$plot_line);
      push(@fv,@fields);
      if( ($#fv+1) >= ($nbn*$nfv) ) { $fvflag = 0; }
    }
  } 
#
#  Write Surfer output file
#
  if( $plot_package =~ /^surfer\b/i ) {
#
      print OUT "\"X, $xunit\" \"Y, $yunit\" \"Z, $zunit\" \"Area, $ftaunit\"";
      for( $i = 0; $i <= $#fvname; $i++ ) {
        if( $fvunit[$i] ) {
          print OUT " \"$fvname[$i], $fvunit[$i]\"";
        } else {
          print OUT " \"$fvname[$i]\"";
        }
      }
      print OUT "\n";
      for( $i = 0; $i < $nbn; $i++ ) {
        $j1 = $i*$nvert;
        $j2 = $j1+1;
        $j3 = $j1+2;
        $j4 = $j1+3;
        $j5 = $j1+4;
        $j6 = $j1+5;
        $j7 = $j1+6;
        $j8 = $j1+7;
        $xs = ($xp[$j1]+$xp[$j2]+$xp[$j3]+$xp[$j4])/8.;
        $xs += ($xp[$j5]+$xp[$j6]+$xp[$j7]+$xp[$j8])/8.;
        $ys = ($yp[$j1]+$yp[$j2]+$yp[$j3]+$yp[$j4])/8.;
        $ys += ($yp[$j5]+$yp[$j6]+$yp[$j7]+$yp[$j8])/8.;
        $zs = ($zp[$j1]+$zp[$j2]+$zp[$j3]+$zp[$j4])/8.;
        $zs += ($zp[$j5]+$zp[$j6]+$zp[$j7]+$zp[$j8])/8.;
        print OUT "$xs $ys $zs $nv[$i]";
        for( $j = 0; $j < $nfv; $j++ ) {
          $k = $j*$nbn + $i;
          print OUT " $fv[$k]";
        }
        print OUT "\n";
      }
  }
#
#  Write Tecplot output file
#
  if( $plot_package =~ /^tecplot\b/i || $plot_package =~ /^tecplot10\b/i ) {
    if( $plot_line == $plot_array[0] ) {
      if( $nf == 0 ) {
        print OUT "TITLE = \"$title\"\n";
        if( $xsurf*$ysurf*$zsurf ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Y, $yunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $xvert == 1 && $yvert == 1 && $zvert == 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Y, $yunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $ifld > 1 && $jfld > 1 && $kfld > 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Y, $yunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $ifld > 1 && $jfld > 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Y, $yunit\" \"Volume, $nvunit\"";
        } elsif( $jfld > 1 && $kfld > 1 ) {
          print OUT "VARIABLES = \"Y, $yunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $kfld > 1 && $ifld > 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $ifld > 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Y, $yunit\" \"Volume, $nvunit\"";
        } elsif( $jfld > 1 ) {
          print OUT "VARIABLES = \"Y, $yunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } elsif( $kfld > 1 ) {
          print OUT "VARIABLES = \"X, $xunit\" \"Z, $zunit\" \"Volume, $nvunit\"";
        } else {
          die "Error: Single Node Plot\n";
        }
        for( $i = 0; $i <= $#fvname; $i++ ) {
          if( $fvunit[$i] ) {
            print OUT "\" $fvname[$i], $fvunit[$i]\"";
          } else {
            print OUT "\" $fvname[$i]\"";
          }
        }
        print OUT "\n";
      }
    }
    if( $ts_opt == 1 ) {
      if( $ta_opt == 1 ) {
        $tsv[$nf] = $tsv[$nf] + $tadj;
      }
      print OUT "ZONE T = \"$tsv[$nf], $tsunit\", STRANDID = 1, SOLUTIONTIME = $tsv[$nf] ";
    } else {
      print OUT "ZONE T = \"$plot_file[$nf]\", ";
    }

    if( $nbn >= 1 ) {
      $nnodes = $nvert*$nbn;
      if( $plot_package =~ /^tecplot10\b/i ) {
        print OUT "N = $nnodes, E = $nbn, DATAPACKING = BLOCK, ZONETYPE = FEBRICK";
      } else {
      print OUT "NODES = $nnodes, ELEMENTS = $nbn, DATAPACKING = BLOCK, ZONETYPE = FEBRICK";
      }  
      if( $nf > 0 ) {
        print OUT " VARSHARELIST = ([1,2,3,4]=1), CONNECTIVITYSHAREZONE = 1\n";
      } else {
        print OUT "\n";
      }
      print OUT "VARLOCATION=([4";
      for( $j = 0; $j < $nfv; $j++ ) {
        $k = $j + 5;
        print OUT ",$k";
      }
      print OUT "]=CELLCENTERED";
      print OUT ")\n";   
    }     
#
#   Tecplot XYZ FETRIANGLE with vertices
#
    if( $nf == 0 ) {
      for( $n = 0; $n < $nbn; $n++ ) {
        for( $m = 0; $m < $nvert; $m++ ) {
          print OUT "$xp[$nc] ";
          $nc++;
        }
        if( $nc > 0 ) {
          print OUT "\n";
        }
      }
      $nc = 0;
      for( $n = 0; $n < $nbn; $n++ ) {
        for( $m = 0; $m < $nvert; $m++ ) {
          print OUT "$yp[$nc] ";
          $nc++;
        }
        if( $nc > 0 ) {
          print OUT "\n";
        }
      }
      $nc = 0;
      for( $n = 0; $n < $nbn; $n++ ) {
        for( $m = 0; $m < $nvert; $m++ ) {
          print OUT "$zp[$nc] ";
          $nc++;
        }
        if( $nc > 0 ) {
          print OUT "\n";
        }
      }
#
      $nc = 0;
      for( $i = 0; $i < $nbn; $i++ ) {
        $nc++;
        print OUT "$nv[$i] ";
        if( $nc == 10 ) {
          print OUT "\n";
          $nc = 0;
        }
      }
      if( $nc > 0 ) {
          print OUT "\n";
      } 
    }          
#
      for( $n = 0; $n < $nfv; $n++ ) {
        $nc = 0;
        for( $i = 0; $i < $nbn; $i++ ) {
          $nc++;
          $k = $n*$nbn + $i;
          print OUT "$fv[$k] ";
          if( $nc == 10 ) {
            print OUT "\n";
            $nc = 0;
          }
        }
        if( $nc > 0 ) {
        print OUT "\n";
        }
    }
#
#   Read Tecplot connectivity list from the 'connect_bh' file
#
    if( $nf == 0 ) {
      open( CONNECT,"connect_bh" ) || die "Error: Unable to Tecplot Fracture Connectivity List File: connect_bh.\n";
      @connect_array = <CONNECT>;
#
#     Loop over lines in connect file
#
      foreach $connect_line (@connect_array) {
#
#       Remove return from line
#
        chomp( $connect_line );
#
#      Remove leading blank spaces from line
#
        $connect_line =~ s/^\s+//;
        print OUT "$connect_line\n";
      }
      close( CONNECT );
    }
  }
#
#  Write VTK output file
#
  $nnodes = $nvert*$nbn;
  if( $plot_package =~ /^vtk\b/i ) {
    if( $plot_line == $plot_array[0] ) {
      print OUT "# vtk DataFile Version 2.0\n";
      print OUT "Unstructured Grid\n";
      print OUT "ASCII\n";
      print OUT "DATASET UNSTRUCTURED_GRID\n";
      print OUT "POINTS $nnodes float\n";
    }
#
#   VTK Unstructured Grid Triangle
#
    $nc = 0;
    for( $n = 0; $n < $nbn; $n++ ) {
      for( $m = 0; $m < $nvert; $m++ ) {
        print OUT "$xp[$nc] $yp[$nc] $zp[$nc] ";
        $nc++;
      }
      if( $nc > 0 ) {
        print OUT "\n";
      }
    }
#
#   Read VTK connectivity list from the 'connect' file and write out cells
#
    open( CONNECT,"connect_bh" ) || die "Error: Unable to Tecplot Connectivity List File: connect_bh.\n";
    $nsize = $nvert*9;
    print OUT "\n";
    print OUT "CELLS $nvert $nsize\n";
    @connect_array = <CONNECT>;
#
#   Loop over lines in connect file
#
    foreach $connect_line (@connect_array) {
#
#     Remove return from line
#
      chomp( $connect_line );
#
#    Remove leading blank spaces from line
#
      $connect_line =~ s/^\s+//;
      @connect_points = split(/\s+/,$connect_line);
      print OUT "$nvert ";
      for( $i = 0; $i <= $#connect_points; $i++ ) {
        $c_one = $connect_points[$i]-1;
        print OUT "$c_one ";
      }
      print OUT "\n"; 
    }
    close( CONNECT );
# 
# Write out cell types (VTK_Hexahedron=12 for 3D)
#
    print OUT "\n";
    print OUT "CELL_TYPES $nvert\n";    
#    if( $ifld > 1 && $jfld > 1 && $kfld > 1 ) {
      for( $m = 0; $m < $nvert; $m++ ) {
        print OUT "12\n";
      }
#    } elsif( ($ifld > 1 && $jfld > 1) || ($jfld > 1 && $kfld > 1) || ($kfld > 1 && $ifld > 1)  ) {
#      for( $m = 0; $m < $nvert; $m++ ) {
#        print OUT "9\n";
#      }
#    } else {
#      for( $m = 0; $m < $nvert; $m++ ) {
#        print OUT "12\n";
#      }
#    }
# 
# Write out field data as cell_data (cell centered)
#
    print OUT "\n";
    print OUT "CELL_DATA $nvert\n";
    print OUT "SCALARS Cell_Volume float\n";
    print OUT "LOOKUP_TABLE default\n";
    $nc = 0;
    for( $i = 0; $i < $nbn; $i++ ) {
      if( $ixp[$i] != 0 ) {
        $nc++;
        print OUT "$nv[$i] ";
        if( $nc == 10 ) {
          print OUT "\n";
          $nc = 0;
        }
      }
    }
    if( $nc > 0 ) {
      print OUT "\n";
    }
    print OUT "\n";
    for( $n = 0; $n <= $#fvname; $n++ ) {
      if( $fvunit[$n] ) {
#        @vtk_name = split(/\s+/,$fvname[$n]);
        $fvname[$n] =~ s/ /_/g;
        print OUT "SCALARS $fvname[$n]_$fvunit[$n] float\n";
        print OUT "LOOKUP_TABLE default\n";
      }else { 
        $fvname[$n] =~ s/ /_/g;
        print OUT "SCALARS $fvname[$n] float\n";
        print OUT "LOOKUP_TABLE default\n";
      }
      $nc = 0;
      for( $i = 0; $i < $nbn; $i++ ) {
        if( $ixp[$i] != 0 ) {
          $nc++;
          $k = $n*$nbn + $i;
          print OUT "$fv[$k] ";
          if( $nc == 10 ) {
            print OUT "\n";
            $nc = 0;
          }
        }
      }
      if( $nc > 0 ) {
        print OUT "\n";
      }
      print OUT "\n";
    }      
  }
  close( PLOT );
  @xp = ();
  @yp = ();
  @zp = ();
  @nv = ();
  @fv = ();
  if( $nf < $#plot_file ) {
    @fvname = ();
    @fvunit = ();
  }
}
close( OUT );
