#!/usr/bin/perl
# subDomain.pl 1322 2020-03-17 21:06:11Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
#
#  Print banner.
#
print("\nWelcome to subDomain ...\n\n");
print("This perl program transforms a STOMP zonation file into\n");
print("a subdomain zonation file.\n");
print("Entries not made on the command line will be prompted.\n\n");
#
#  Prompt user for domain dimensions.
#
print("Domain Dimensions?\n");
$domain = <STDIN>;
chomp( $domain );
@fields = split(/\s+|,/,$domain);
$ifld = $fields[0];
$jfld = $fields[1];
$kfld = $fields[2];
#print "IFLD = $ifld, JFLD = $jfld, KFLD = $kfld\n";
if( $ifld < 1 ) { die "\nDomain Error: IFLD < 1\n" };
if( $jfld < 1 ) { die "\nDomain Error: JFLD < 1\n" };
if( $kfld < 1 ) { die "\nDomain Error: KFLD < 1\n" };
print("Lower Subdomain Dimensions?\n");
$domain = <STDIN>;
chomp( $domain );
@fields = split(/\s+|,/,$domain);
$ilsd = $fields[0];
$jlsd = $fields[1];
$klsd = $fields[2];
#print "IL = $ilsd, JL = $jlsd, KL = $klsd\n";
if( $ilsd < 1 ) { die "\nLower Subdomain Error: IL < 1\n" };
if( $jlsd < 1 ) { die "\nLower Subdomain Error: JL < 1\n" };
if( $klsd < 1 ) { die "\nLower Subdomain Error: KL < 1\n" };
if( $ilsd > $ifld ) { die "\nLower Subdomain Error: IL > IFLD\n" };
if( $jlsd > $jfld ) { die "\nLower Subdomain Error: JL > JFLD\n" };
if( $klsd > $kfld ) { die "\nLower Subdomain Error: KL > KFLD\n" };
print("Upper Subdomain Dimensions?\n");
$domain = <STDIN>;
chomp( $domain );
@fields = split(/\s+|,/,$domain);
$iusd = $fields[0];
$jusd = $fields[1];
$kusd = $fields[2];
#print "IU = $iusd, JU = $jusd, KU = $kusd\n";
if( $iusd < $ilsd ) { die "\nUpper Subdomain Error: IU < IL\n" };
if( $jusd < $jlsd ) { die "\nUpper Subdomain Error: JU < JL\n" };
if( $kusd < $klsd ) { die "\nUpper Subdomain Error: KU < KL\n" };
if( $iusd > $ifld ) { die "\nUpper Subdomain Error: IU > IFLD\n" };
if( $jusd > $jfld ) { die "\nUpper Subdomain Error: JU > JFLD\n" };
if( $kusd > $kfld ) { die "\nUpper Subdomain Error: KU > KFLD\n" };
#
#  Domain file name as first argument or prompt user.
#
if( $ARGV[0] ) {
  $domain_file = $ARGV[0];
  chomp( $domain_file );
  open( DOMAIN,$domain_file) || die "Error: Unable to Open Domain File: $domain_file.\n";
  if( $ARGV[1] ) {
    $subdomain_file = $ARGV[1];
    chomp( $subdomain_file );
    open( SUBDOMAIN,">$subdomain_file") || die "Error: Unable to Open Subdomain File: $subdomain_file.\n";
#
#  No second argument; ask user for plotting-package input file name.
#
  } else {
    $stops = 0;
    do {
      print("Subdomain File Name?\n");
      $subdomain_file = <STDIN>;
      chomp( $subdomain_file );
      if( open( SUBDOMAIN,">$subdomain_file" ) ) {
        $stops = 1;
      } else {
        print("Error: Unable to Open Subdomain File: $subdomain_file. Try again!\n\n");
      }
    } until $stops != 0;
  }
} else {
#
#  No first or second argument; ask user for domain and subdomain file names.
#
  $stops = 1;
  do {
    print("Domain File Name?\n");
    $domain_file = <STDIN>;
    chomp( $domain_file );
    if( -r $domain_file ) {
      $stops = 1;
    } elsif( -B $in_file )  {
      $stops = 0;
      print("Error: Domain File is Binary: $domain_file.\n");
    } else {
      $stops = 0;
      print("Error: Domain File is Unreadable: $domain_file.\n");
    }
    if( $stops == 0 ) { print("Try again!\n\n"); }
  } until $stops != 0;
  open( DOMAIN,$domain_file) || die "Error: Unable to Open Domain File: $domain_file.\n";
  $stops = 0;
  do {
    print("Subdomain File Name?\n");
    $subdomain_file = <STDIN>;
    chomp( $subdomain_file );
    if( open( SUBDOMAIN,">$subdomain_file" ) ) {
      $stops = 1;
    } else {
      print("Error: Unable to Open Subdomain File: $subdomain_file. Try again!\n\n");
    }
  } until $stops != 0;
}
#
#  Initialize
#
$nfld = 0;
#
#  Read domain file into domain array
#
@domain_array = <DOMAIN>;
#
#  Loop over lines in domain file
#
foreach $domain_line (@domain_array) {
#
#  Remove return from line
#
  chomp( $domain_line );
#
#  Remove leading blank spaces from line
#
  $domain_line =~ s/^\s+//;
  @fields = split(/\s+|,/,$domain_line);
  push( @iz,@fields );
  $nfld += $#fields + 1;
}
$ijkfld = $ifld*$jfld*$kfld;
$ijfld = $ifld*$jfld;
if( $nfld != $ijkfld ) { die "\nDomain Error: NFLD ($nfld) != IFLD*JFLD*KFLD ($ijkfld)\n" };
#
#  Load subdomain array
#
for( $k = $klsd; $k <= $kusd; $k++ ) {
  for( $j = $jlsd; $j <= $jusd; $j++ ) {
    for( $i = $ilsd; $i <= $iusd; $i++ ) {
      $n = ($k-1)*$ijfld + ($j-1)*$ifld + $i;
      $n--;
      push( @subdomain_array,$iz[$n] );
    }
  }
}
#
#  Write subdomain
#
for( $n = 0; $n <= $#subdomain_array; $n++ ) {
  push( @line,$subdomain_array[$n] );
  if( $#line >= 9 ) {
    print SUBDOMAIN "@line\n";
    @line = ();
  }
}
print SUBDOMAIN "@line\n";

