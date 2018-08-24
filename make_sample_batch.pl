#!/usr/bin/perl

my $REF_FILE = "/disk01/usr5/matanaka/proton_decay/skdetsim-v13p90_dec17";
my $OUTPUT_FILE = "/disk01/usr5/matanaka/proton_decay/samples";
my $INPUT_FILE = "/disk01/usr5/matanaka/proton_decay/atmpd_17a/src/analysis/ndecay";
{
  my $RUN_NUM_BEGIN=$ARGV[0];
  my $RUN_NUM_END=$ARGV[1];
  my $mode=$ARGV[2];
  my $CARD_NAME=$ARGV[3];
  $m=$RUN_NUM_BEGIN;

  system "echo Running sample $RUN_NUM_BEGIN to $RUN_NUM_END";
  system "echo signal is $mode";

  use File::Path;

  system "mkdir -p $OUTPUT_FILE/$mode";
  system "mkdir -p $OUTPUT_FILE/$mode/script";
  system "mkdir -p $OUTPUT_FILE/$mode/log";
  system "mkdir -p $OUTPUT_FILE/$mode/err";
  system "mkdir -p $OUTPUT_FILE/$mode/file";

  while ($m <=$RUN_NUM_END){
    my $input_skdetsim = sprintf("$INPUT_FILE/$mode/output/vector%d",$m);
    my $output_skdetsim=sprintf("../samples/$mode/file/output_skdetsim.%d",$m);
    my $input_apfit = $output_skdetsim;
    my $output_apfit=sprintf("../samples/$mode/file/output_apfit.%d",$m);
    my $Script=sprintf("$OUTPUT_FILE/$mode/script/script.%d.csh",$m);
    my $LogFile=sprintf("$OUTPUT_FILE/$mode/log/logfile.%d.log",$m);
    my $ErrFile=sprintf("$OUTPUT_FILE/$mode/err/errfile.%d.log",$m);
    WriteScriptFile($Script,$m,$input_skdetsim,$output_skdetsim,$input_apfit,$output_apfit,$CARD_NAME);
    $cmd="qsub -q atmpd -o $LogFile -e $ErrFile $Script";
    system "echo $cmd";
    system $cmd;
    $m++;
  }
}

#------------------------#                                                                                                       
sub WriteScriptFile() {
#------------------------#                                                                                                       
  my $file   = $_[0];
  my $run = $_[1];
  my $in_skdetsim = $_[2];
  my $out_skdetsim = $_[3];
  my $in_apfit = $_[4];
  my $out_apfit = $_[5];
  my $card_name = $_[6];
  open (SCRIPT,">$file");
  print SCRIPT "#!/bin/csh -f\n";
  print SCRIPT "cd $REF_FILE\n";
  print SCRIPT "source setup.sh\n";
  print SCRIPT "date\n";
  print SCRIPT "hostname\n";
  print SCRIPT "./skdetsim.sh $card_name $out_skdetsim $in_skdetsim\n"; 
  print SCRIPT "./apfit_fc.sh $in_apfit $out_apfit\n";
  close SCRIPT;
  $cmdchmod = "chmod 755 $file";
  system $cmdchmod;
}
