#!/usr/bin/perl

# Please see the LICENSE accompanying this distribution for 
# details.  Report all bugs to nhomer@cs.ucla.edu or 
# srma-help@lists.sourceforge.net.  For documentation, use the
# -man option.

# TODO:
# - guarantee tmp files do not collide
# - documentation
# - documentation: split size == 0 means it wont split

use strict;
use warnings FATAL => qw( all );
use File::Path;
use XML::Simple; 
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Path; # for directory creation
use Cwd;

my %QUEUETYPES = ("SGE" => 0, "PBS" => 1);
my %SPACE = ("NT" => 0, "CS" => 1);
my $FAKEQSUBID = 0;

use constant {
	OPTIONAL => 0,
	REQUIRED => 1,
	BREAKLINE => "************************************************************\n",
	MERGE_LOG_BASE => 8,
	QSUBNOJOB => "QSUBNOJOB"
};

my $config;
my ($man, $print_schema, $help, $quiet, $dryrun) = (0, 0, 0, 0, 0);
my $version = "0.1.1";

GetOptions('help|?' => \$help, 
	man => \$man, 
	'schema' => \$print_schema, 
	'quiet' => \$quiet, 
	'dryrun' => \$dryrun,
	'config=s' => \$config)
	or pod2usage(1);
Schema() if ($print_schema);
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if ($help or !defined($config));

if(!$quiet) {
	print STDOUT BREAKLINE;
}

# Read in from the config file
my $xml = new XML::Simple;
my $data = $xml->XMLin("$config");

# Validate data
ValidateData($data);

# Submit jobs
CreateJobs($data, $quiet, $dryrun);

if(!$quiet) {
	print STDOUT BREAKLINE;
}

sub Schema {
	# print the schema
	my $schema = <<END;
	<?xml version="1.0"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:element name="srmaConfig">
	<xs:complexType>
	  <xs:sequence>
		<xs:element name="srmaOptions">
		  <xs:complexType>
			<xs:sequence>
			  <xs:element name="srmaBin" type="directoryPath"/>
			  <xs:element name="javaBin" type="directoryPath"/>
			  <xs:element name="qsubBin" type="directoryPath"/>
			  <xs:element name="referenceFasta" type="filePath" use="required/>
				<xs:complexType>
				  <xs:attribute name="splitSize" type="positiveInteger" use="required"/>
				</xs:complexType>
			  </xs:element>
			  <xs:element name="offset" type="xs:integer" use="required"/>
			  <xs:element name="minimumAlleleFrequency" type="xs:double" use="required"/>
			  <xs:element name="minimumAlleleCoverage" type="xs:integer" use="required"/>
			  <xs:element name="inputBAMFile" type="filePath" use="required">
			  <xs:element name="outputBAMFile" type="filePath" use="required">
			  <xs:element name="runDirectory" type="directoryPath" use="required">
			  <xs:element name="tmpDirectory" type="directoryPath" use="required"/>
			  <xs:element name="queueType" use="required">
				<xs:simpleType>
				  <xs:restriction base="xs:string">
					<xs:enumeration value="SGE"/>
					<xs:enumeration value="PBS"/>
				  </xs:restriction>
				</xs:simpleType>
			  </xs:element>
			  <xs:element name="cleanUpTmpDirectory" type="xs:integer"/>
			  <xs:element name="javaArgs" type="xs:string"/>
			  <xs:element name="qsubArgs" type="xs:string"/>
			</xs:sequence>
		  </xs:complexType>
		</xs:element>
		<xs:element name="samOptions">
		  <xs:complexType>
			<xs:sequence>
			  <xs:element name="picardBin" type="directoryPath"/>
			  <xs:element name="javaArgs" type="xs:string"/>
			  <xs:element name="qsubArgs" type="xs:string"/>
			</xs:sequence>
		  </xs:complexType>
		</xs:element>
	  </xs:sequence>
	</xs:complexType>
  </xs:element>
  <xs:simpleType name="filePath">
	<xs:restriction base="xs:string">
	  <xs:pattern value="\\S+"/>
	</xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="directoryPath">
	<xs:restriction base="xs:string">
	  <xs:pattern value="\\S+/"/>
	</xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="positiveInteger">
	<xs:restriction base="xs:integer">
	  <xs:minInclusive value="1"/>
	</xs:restriction>
  </xs:simpleType>
</xs:schema>
END

	print STDOUT $schema;
	exit 1;
}

sub ValidateData {
	my $data = shift;

	# global options
	die("The global options were not found.\n") unless (defined($data->{'srmaOptions'})); 
	ValidatePath($data->{'srmaOptions'},         'srmaBin',                                  OPTIONAL); 
	ValidatePath($data->{'srmaOptions'},         'javaBin',                                  OPTIONAL); 
	ValidateFile($data->{'srmaOptions'},         'javaArgs',							     OPTIONAL);
	ValidatePath($data->{'srmaOptions'},         'qsubBin',                                  OPTIONAL); 
	ValidateOptions($data->{'srmaOptions'},      'queueType',          \%QUEUETYPES,         REQUIRED);
	ValidateOption($data->{'srmaOptions'},	     'offset',									 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'minimumAlleleFrequency',					 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'minimumAlleleCoverage',					 OPTIONAL);
	ValidateFile($data->{'srmaOptions'},         'referenceFasta',                           REQUIRED);
	ValidatePath($data->{'srmaOptions'},         'runDirectory',                             REQUIRED); 
	ValidatePath($data->{'srmaOptions'},         'tmpDirectory',                             REQUIRED); 
	ValidateFile($data->{'srmaOptions'},         'inputBAMFile',							 REQUIRED);
	ValidateFile($data->{'srmaOptions'},         'outputBAMFile',							 REQUIRED);

	die "Attribute splitSize required with referenceFasta.\n" if (!defined($data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'}));
	die "Attribute splitSize must be greater than or equal tozero.\n" if ($data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'} < 0);

	# picard
	if(defined($data->{'samOptions'})) {
		ValidatePath($data->{'samOptions'},       'picardBin',                                OPTIONAL); 
		ValidateOption($data->{'samOptions'},     'cleanUpTmpDirectory',                      OPTIONAL);
		ValidateOption($data->{'samOptions'},     'qsubArgs',                                 OPTIONAL);
	}
}

sub ValidateOption {
	my ($hash, $option, $required) = @_;

	return 1 if (defined($hash->{$option}));
	return 0 if (OPTIONAL == $required);

	die("Option '$option' was not found.\n");
}

sub ValidatePath {
	my ($hash, $option, $required) = @_; 

	if(0 != ValidateOption($hash, $option, $required) and $hash->{$option} !~ m/\S+\//) { # very liberal check
		die("Option '$option' did not give a valid path.\n");
	}
}

sub ValidateFile {
	my ($hash, $option, $required) = @_;

	if(0 != ValidateOption($hash, $option, $required) and $hash->{$option} !~ m/\S+/) { # very liberal check
		die("Option '$option' did not give a valid file name.\n");
	}
}

sub ValidateOptions {
	my ($hash, $option, $values, $required) = @_;

	if(0 != ValidateOption($hash, $option, $required) and !defined($values->{$hash->{$option}})) {
		die("The value '".($hash->{$option})."' for option '$option' was not valid.\n");
	}
}

sub getGenomeInfo {
	my ($fasta, $genomeInfo) = @_;

	open(FH, "$fasta.fai") || die("Could not open FASTA index: $fasta.fai");
	while(defined(my $line = <FH>)) {
		my @a = split(/\s+/, $line);
		$genomeInfo->{$a[0]} = $a[1];
	}
	close(FH);

	# TODO
	return 0;
}

sub CreateJobs {
	my ($data, $quiet, $dryrun) = @_;

	my @srmaJobIDs = ();
	my @srmaOutputIDs = ();
	my @samJobIDs = ();
	my @samOutputIDs = ();

	# Create directories - Error checking...
	mkpath([$data->{'srmaOptions'}->{'runDirectory'}],    ($quiet) ? 0 : 1, 0755);
	mkpath([$data->{'srmaOptions'}->{'tmpDirectory'}],    ($quiet) ? 0 : 1, 0755);

	CreateJobsSRMA($data, $quiet, $dryrun, \@srmaJobIDs, \@srmaOutputIDs);
	if(0 < scalar(@srmaJobIDs)) {
		CreateJobsSAM($data, $quiet, $dryrun, \@srmaJobIDs, \@srmaOutputIDs);
	}
}

sub CreateRunFile {
	my ($data, $type, $outputID) = @_;

	if(0 < length($outputID)) {
		return $data->{'srmaOptions'}->{'runDirectory'}."$type.$outputID.sh";
	}
	else {
		return $data->{'srmaOptions'}->{'runDirectory'}."$type.sh";
	}
}

sub CreateTmpOutputFile {
	my ($data, $type, $outputID) = @_;
	return $data->{'srmaOptions'}->{'tmpDirectory'}."$type.$outputID.bam";
}

sub CreateJobsSRMA {
	my ($data, $quiet, $dryrun, $qsubIDs, $outputIDs) = @_;
	my @read_files = ();

	if(0 == $data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'}) { # do not split
		my $runFile = CreateRunFile($data, 'srma', "");
		my $cmd = "";
		$cmd = $data->{'srmaOptions'}->{'javaBin'}."java";
		if(defined($data->{'srmaOptions'}->{'javaArgs'})) {
			$cmd .= " ".$data->{'srmaOptions'}->{'javaArgs'};
			if($data->{'srmaOptions'}->{'javaArgs'} !~ m/-Xmx/) {
				$cmd .= " -Xmx2g";
			}
		}
		else {
			$cmd .= " -Xmx2g";
		}
		$cmd .= " -jar ".$data->{'srmaOptions'}->{'srmaBin'}."srma.jar";
		$cmd .= " I=".$data->{'srmaOptions'}->{'inputBAMFile'};
		$cmd .= " O=".$data->{'srmaOptions'}->{'outputBAMFile'};
		$cmd .= " R=".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
		$cmd .= " OFFSET=".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
		$cmd .= " MINIMUM_ALLELE_FREQUENCY=".$data->{'srmaOptions'}->{'minimumAlleleFrequency'} if(defined($data->{'srmaOptions'}->{'minimumAlleleFrequency'}));
		$cmd .= " MINIMUM_ALLELE_COVERAGE=".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
		$cmd .= " QUIET=true";

		# Submit the job
		my @a = (); # empty array for job dependencies
		my $qsubID = SubmitJob($runFile, $quiet, 0, $dryrun, $cmd, $data, 'srmaOptions', 'srma', \@a);
		# No need to merge
	}
	else {
		my $start = 1;
		my $splitSize = $data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'};
		my %genomeInfo = ();

		getGenomeInfo($data->{'srmaOptions'}->{'referenceFasta'}->{'content'}, \%genomeInfo);

		foreach my $chrName (keys %genomeInfo) {
			my $chrSize = $genomeInfo{$chrName};
			for(my $start=1;$start <= $chrSize;$start+=$splitSize) {
				my $end = $start + $splitSize - 1;
				if($chrSize < $end) { $end = $chrSize; }
				my $outputID = "$chrName\_$start\-$end";
				my $runFile = CreateRunFile($data, 'srma', $outputID);
				my $outputFile = CreateTmpOutputFile($data, 'srma', $outputID);

				my $cmd = "";
				$cmd = $data->{'srmaOptions'}->{'javaBin'}."java";
				if(defined($data->{'samOptions'}->{'javaArgs'})) {
					$cmd .= " ".$data->{'samOptions'}->{'javaArgs'};
					if($data->{'samOptions'}->{'javaArgs'} !~ m/-Xmx/) {
						$cmd .= " -Xmx2g";
					}
				}
				else {
					$cmd .= " -Xmx2g";
				}
				$cmd .= " -jar ".$data->{'srmaOptions'}->{'srmaBin'}."srma.jar";
				$cmd .= " I=".$data->{'srmaOptions'}->{'inputBAMFile'};
				$cmd .= " O=$outputFile";
				$cmd .= " R=".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
				$cmd .= " OFFSET=".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
				$cmd .= " MINIMUM_ALLELE_FREQUENCY=".$data->{'srmaOptions'}->{'minimumAlleleFrequency'} if(defined($data->{'srmaOptions'}->{'minimumAlleleFrequency'}));
				$cmd .= " MINIMUM_ALLELE_COVERAGE=".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
				$cmd .= " RANGE=\"$chrName\:$start-$end\"";
				$cmd .= " QUIET=true";

				# Submit the job
				my @a = (); # empty array for job dependencies
				my $qsubID = SubmitJob($runFile, $quiet, 0, $dryrun, $cmd, $data, 'srmaOptions', $outputID, \@a);
				push(@$qsubIDs, $qsubID) if (QSUBNOJOB ne $qsubID);
				push(@$outputIDs, $outputID);
			}
		}
	}
}


sub CreateJobsSAM {
	my ($data, $quiet, $dryrun, $dependentQsubIDs, $dependentOutputIDs) = @_;
	my @qsubIDs = ();
	my ($cmd, $run_file, $outputID, $qsub_id);

	my $type = 'picard';
	my @tmpOutputBAMFiles = ();

	if(!defined($data->{'samOptions'})) {
		# Delete the dependent jobs
		for(my $i=0;$i<scalar($dependentQsubIDs);$i++) {
			my $qsubID = $dependentQsubIDs->[$i];
			`qdel $qsubID`;
		}
		die("samOptions must be defined");
	}

	# Get BAM file names
	for(my $i=0;$i<scalar(@$dependentOutputIDs);$i++) {
		my $outputID = $dependentOutputIDs->[$i];
		my $tmpBAMFile = CreateTmpOutputFile($data, 'srma', $outputID);
		push(@tmpOutputBAMFiles, $tmpBAMFile);
		push(@qsubIDs, $dependentQsubIDs->[$i]);
	}

	# Merge script(s)
	# Note: there could be too many dependencies, so lets just create dummy jobs to "merge" the dependencies
	my $mergeLevel = 0;
	while(MERGE_LOG_BASE < scalar(@qsubIDs)) { # while we must merge
		$mergeLevel++;
		my $ctr = 0;
		my @curIDs = @qsubIDs;
		@qsubIDs = ();
		for(my $i=0;$i<scalar(@curIDs);$i+=MERGE_LOG_BASE) {
			$ctr++;
			# Get the subset of dependent jobs
			my @dependent_jobs = ();
			for(my $j=$i;$j<scalar(@curIDs) && $j<$i+MERGE_LOG_BASE;$j++) {
				push(@dependent_jobs, $curIDs[$j]);
			}
			# Create the command
			$outputID = "dummy.merge.$mergeLevel.$ctr";
			$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$type.".$outputID.".sh";
			$cmd = "echo \"Merging $mergeLevel / $ctr\"\n";
			$qsub_id = SubmitJob($run_file, $quiet, 1, $dryrun, $cmd, $data, 'samOptions', $outputID, \@dependent_jobs);
			if(QSUBNOJOB ne $qsub_id) {
				push(@qsubIDs, $qsub_id);
			}
			else {
				# currently it must be submitted
				die;
			}
		}
	}

	# Run merge or copy
	if(1 < scalar(@qsubIDs)) {
		$outputID = "merge";
		$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$type.".$outputID.".sh";
		if(!defined($data->{'samOptions'}->{'picardBin'})) { die("Picard bin required") };
		$cmd = $data->{'srmaOptions'}->{'javaBin'}."java";
		if(defined($data->{'samOptions'}->{'javaArgs'})) {
			$cmd .= " ".$data->{'samOptions'}->{'javaArgs'};
			if($data->{'samOptions'}->{'javaArgs'} !~ m/-Xmx/) {
				$cmd .= " -Xmx2g";
			}
		}
		else {
			$cmd .= " -Xmx2g";
		}
		$cmd .= " -jar ".$data->{'samOptions'}->{'picardBin'}."MergeSamFiles.jar";
		foreach my $bam (@tmpOutputBAMFiles) {
			$cmd .= " I=$bam";
		}
		$cmd .= " O=".$data->{'srmaOptions'}->{'outputBAMFile'};
		$cmd .= " SO=coordinate";
		$cmd .= " AS=true";
		$cmd .= " TMP_DIR=".$data->{'srmaOptions'}->{'tmpDirectory'};
		$cmd .= " VALIDATION_STRINGENCY=SILENT";
	}
	else {
		$outputID = "copy";
		$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$type.".$outputID.".sh";
		my $bam = $tmpOutputBAMFiles[0];
		$cmd = "cp -v $bam ".$data->{'srmaOptions'}->{'outputBAMFile'};
	}
	$qsub_id = SubmitJob($run_file , $quiet, 1, $dryrun, $cmd, $data, 'samOptions', $outputID, \@qsubIDs);
	if(QSUBNOJOB ne $qsub_id) {
		push(@qsubIDs, $qsub_id);
	}
	else {
		# currently it must be submitted
		die;
	}

	# Clean up
	if(defined($data->{'srmaOptions'}->{'cleanUpTmpDirectory'}) && 1 == defined($data->{'srmaOptions'}->{'cleanUpTmpDirectory'})) {
		my @a = (); push(@a, $qsub_id); # push the merge/copy
		$outputID = "cleanup.tmpdirectory";
		$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$outputID.sh";
		$cmd = "rm -rv ".$data->{'srmaOptions'}->{'tmpDirectory'};
		$qsub_id = SubmitJob($run_file , $quiet, 1, $dryrun, $cmd, $data, 'srmaOptions', $outputID, \@a);
	}
}

sub SubmitJob {
	my ($run_file, $quiet, $should_depend, $dryrun, $command, $data, $type, $outputID, $dependent_jobIDs) = @_;
	if(0 < length($outputID)) {
		$outputID = "$type.$outputID"; $outputID =~ s/Options//g;
	}
	else {
		$outputID = "$type"; $outputID =~ s/Options//g;
	}

	if(!$quiet) {
		print STDERR "[srma submit] RUNFILE=$run_file\n";
	}
	my $output = <<END_OUTPUT;
run ()
{
	echo "running: \$*" 2>&1;
	eval \$*;
	if test \$? != 0 ; then
	echo "error: while running '\$*'";
	exit 100;
	fi
}
END_OUTPUT
	$output .= "\nrun \"hostname\";\n";
	$output .= "run \"$command\";\n";
	$output .= "exit 0;\n";
	open(FH, ">$run_file") or die("Error.  Could not open $run_file for writing!\n");
	print FH "$output";
	close(FH);

	# Create qsub command
	my $qsub = "";
	$qsub .= $data->{'srmaOptions'}->{'qsubBin'} if defined($data->{'srmaOptions'}->{'qsubBin'});
	$qsub .= "qsub";

	if(0 < scalar(@$dependent_jobIDs) && 1 == $should_depend) {
		$qsub .= " -hold_jid ".join(",", @$dependent_jobIDs)         if ("SGE" eq $data->{'srmaOptions'}->{'queueType'});
		$qsub .= " -W depend=afterok:".join(":", @$dependent_jobIDs) if ("PBS" eq $data->{'srmaOptions'}->{'queueType'});
	}
	$qsub .= " ".$data->{$type}->{'qsubArgs'} if defined($data->{$type}->{'qsubArgs'});
	$qsub .= " -N $outputID -o $run_file.out -e $run_file.err $run_file";

	if(1 == $dryrun) {
		$FAKEQSUBID++;
		print STDERR "[srma submit] NAME=$outputID QSUBID=$FAKEQSUBID\n";
		return $FAKEQSUBID;
	}

	# Submit the qsub command
	my $qsub_id=`$qsub`;
	$qsub_id = "$qsub_id";
	chomp($qsub_id);

	# There has to be a better way to get the job ids (?)
	if($qsub_id =~ m/Your job (\d+)/) {
		$qsub_id = $1;
	}
	die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
	if($qsub_id !~ m/^\S+$/) {
		die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
	}

	if(!$quiet) {
		print STDERR "[srma submit] NAME=$outputID QSUBID=$qsub_id\n";
	}

	return $qsub_id;
}

__END__
=head1 SYNOPSIS

srma.submit.pl [options] 

=head1 OPTIONS

=over 8

=item B<-help>
Print a brief help message and exits.

=item B<-schema>
Print the configuration XML schema.

=item B<-man>
Prints the manual page and exits.

=item B<-quiet>
Do not print any submit messages.

=item B<-dryrun>
Do everything but submit the jobs.

=item B<-config>
The XML configuration file.

=back

=head1 DESCRIPTION

B<srma.submit.pl> will create the necessary shell scripts for B<SRMA> to be
run on a supported cluster (SGE and PBS).  It will also submit each script
supporting job dependencies.  The input to B<srma.submit.pl> is an XML
configuration file.  To view the schema for this file, please use the 
I<-schema> option.

Please report all bugs to nhomer@cs.ucla.edu or srma-help@lists.sourceforge.net.

=head1 REQUIREMENTS

This script requires the XML::Simple library, which can be found at:
http://search.cpan.org/dist/XML-Simple

=cut
