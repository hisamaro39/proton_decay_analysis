#!/usr/bin/perl

my $REF_FILE = "/disk02/usr6/matanaka/proton_decay/proton_decay_analysis";
my $OUTPUT_FILE = "/disk02/usr6/matanaka/proton_decay/proton_decay_analysis/output_batch";
{

  my $refArray = @ARGV;

  if($refArray!=4) {
    print "./run_batch.pl [Separation] [INPUT] [MODE] [skx]\n";
    exit;
  }
  my $RUN_NUM_BEGIN=0;
  my $RUN_NUM_END=$ARGV[0];
  my $INPUT=$ARGV[1];
  my $MODE=$ARGV[2];
  my $period=$ARGV[3];
  my $CARD_NAME="./18a_${period}.card";
  $m=$RUN_NUM_BEGIN;

  print "Input is ", $INPUT, "\n";
  print "Mode is ", $MODE, "\n";
  
  my $DIR = $INPUT + '_' + $MODE;

  use File::Path;

  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/script";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/log";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/err";
  system "mkdir -p $OUTPUT_FILE/${INPUT}_$MODE/${period}/file";

  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/script/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/log/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/err/*";
  system "rm $OUTPUT_FILE/${INPUT}_$MODE/${period}/file/*";

  while ($m <$RUN_NUM_END){
    my $Script=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/script/script.%d.csh",$m);
    my $LogFile=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/log/logfile.%d.log",$m);
    my $ErrFile=sprintf("$OUTPUT_FILE/${INPUT}_$MODE/${period}/err/errfile.%d.log",$m);
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
  print SCRIPT "./build_osc_ntuple $per $in $mode $card $run $blossom 1\n"; 
  close SCRIPT;
  $cmdchmod = "chmod 755 $file";
  system $cmdchmod;
}
