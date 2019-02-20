#!/usr/bin/perl

my $REF_FILE = "/disk02/usr6/matanaka/proton_decay/proton_decay_analysis";
my $OUTPUT_FILE = "/disk02/usr6/matanaka/proton_decay/proton_decay_analysis/output_batch";
{

  my $refArray = @ARGV;

  if($refArray!=3) {
    print "./run_batch.pl [INPUT] [MODE] [skx]\n";
    exit;
  }
  my $RUN_NUM_BEGIN=2500;
  my $RUN_NUM_END=3000;
  my $INPUT=$ARGV[0];
  my $MODE=$ARGV[1];
  my $period=$ARGV[2];
  $m=$RUN_NUM_BEGIN;

  print "Input is ", $INPUT, "\n";
  print "Mode is ", $MODE, "\n";
  
  my $DIR = $INPUT + '_' + $MODE;

  use File::Path;

  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/script_one_by_one";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/log_one_by_one";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/err_one_by_one";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/file_one_by_one";

  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/script_one_by_one/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/log_one_by_one/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/err_one_by_one/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/file_one_by_one/*";

  while ($m <$RUN_NUM_END){
    my $Script=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/script_one_by_one/script.%d.csh",$m);
    my $LogFile=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/log_one_by_one/logfile.%d.log",$m);
    my $ErrFile=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/err_one_by_one/errfile.%d.log",$m);
    my $CARD_NAME="./output_batch/Card/18a_${period}.${m}.card";
    WriteScriptFile($Script,$m,$INPUT,$MODE,$CARD_NAME,$RUN_NUM_END,$period);
    $cmd="qsub -q atmpd -o $LogFile -e $ErrFile $Script";
    system $cmd;
    $m++;
  }
}

#------------------------#                                                                                                       
sub WriteScriptFile() {
#------------------------#                                                                                                       
  my $file   = $_[0];
  my $run = $_[1];
  my $in = $_[2];
  my $mode = $_[3];
  my $card = $_[4];
  my $blossom = $_[5];
  my $per = $_[6];
  open (SCRIPT,">$file");
  print SCRIPT "#!/bin/csh -f\n";
  print SCRIPT "/bin/csh \n";
  print SCRIPT "cd $REF_FILE\n";
  print SCRIPT "source setup.sh\n";
  print SCRIPT "date\n";
  print SCRIPT "hostname\n";
  print SCRIPT "./build_osc_ntuple $per $in $mode $card 0 1 1\n"; 
  close SCRIPT;
  $cmdchmod = "chmod 755 $file";
  system $cmdchmod;
}
