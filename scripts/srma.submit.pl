#!/usr/bin/perl

# Please see the LICENSE accompanying this distribution for 
# details.  Report all bugs to nhomer@cs.ucla.edu or 
# srma-help@lists.sourceforge.net.  For documentation, use the
# -man option.

# TODO:
# - guarantee tmp files do not collide
# - documentation
# 	- documentation: split size == 0 means it wont split
# - could split not by regions but by the # of records using the ".bam.bai" file (?)
# - check for presence of ref.fa.dict and in.bam.bai

use strict;
use warnings FATAL => qw( all );
use File::Path;
use XML::Simple; 
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Path; # for directory creation
use Cwd;
use POSIX qw/ceil/; # for ceil command

my %QUEUETYPES = ("SGE" => 0, "PBS" => 1);
my %SPACE = ("NT" => 0, "CS" => 1);
my %STARTSTEP = ("srma" => 0, "sam" => 1);
my $FAKEQSUBID = 0;

use constant {
	OPTIONAL => 0,
	REQUIRED => 1,
	BREAKLINE => "************************************************************\n",
	MERGE_LOG_BASE => 16,
	QSUBNOJOB => "QSUBNOJOB"
};

my $config;
my ($man, $print_schema, $help, $quiet, $start_step, $dryrun) = (0, 0, 0, 0, "srma", 0);
my $version = "0.1.7";

GetOptions('help|?' => \$help, 
	man => \$man, 
	'schema' => \$print_schema, 
	'quiet' => \$quiet, 
	'startstep=s' => \$start_step,
	'dryrun' => \$dryrun,
	'config=s' => \$config)
	or pod2usage(1);
Schema() if ($print_schema);
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(1) if ($help or !defined($config));

if(!defined($STARTSTEP{$start_step})) {
	print STDERR "Error. Illegal value to the option \"$start_step\" -startstep.\n";
	pod2usage(1);
}
$start_step = $STARTSTEP{$start_step};

if(!$quiet) {
	print STDOUT BREAKLINE;
}

# Read in from the config file
my $xml = new XML::Simple;
my $data = $xml->XMLin("$config");
my @qsubGlobalIDs = ();

# Validate data
ValidateData($data);

# Submit jobs
CreateJobs($data, $quiet, $start_step, $dryrun);

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
			  <xs:element name="srmaJar" type="directoryPath"/>
			  <xs:element name="srmaBinary" type="directoryPath"/>
			  <xs:element name="javaBin" type="directoryPath"/>
			  <xs:element name="qsubBin" type="directoryPath"/>
			  <xs:element name="referenceFasta" type="filePath" use="required"/>
				<xs:complexType>
				  <xs:attribute name="splitSize" type="positiveInteger" use="optional"/>
				</xs:complexType>
			  </xs:element>
			  <xs:element name="offset" type="xs:integer" use="optional"/>
			  <xs:element name="minimumAlleleProbability" type="xs:double" use="optional"/>
			  <xs:element name="minimumAlleleCoverage" type="xs:integer" use="optional"/>
			  <xs:element name="maximumTotalCoverage" type="xs:integer" use="optional"/>
			  <xs:element name="inputBAMFile" type="filePath" use="required"/>
			  <xs:element name="outputBAMFile" type="filePath" use="required"/>
			  <xs:element name="runDirectory" type="directoryPath" use="required"/>
			  <xs:element name="tmpDirectory" type="directoryPath" use="required"/>
			  <xs:element name="queueType" type="xs:string "use="required">
				<xs:complexType>
				  <xs:attribute name="arrayJob" type="xs:boolean" default="false" use="optional"/>
				</xs:complexType>
			  </xs:element>
			  <xs:element name="cleanUpTmpDirectory" type="xs:integer"/>
			  <xs:element name="javaArgs" type="xs:string"/>
			  <xs:element name="qsubArgs" type="xs:string"/>
			  <xs:element name="range" type="xs:string"/>
			  <xs:element name="ranges" type="filePath"/>
			  <xs:element name="correctBases" type="xs:string"/>
			  <xs:element name="numThreads" type="xs:integer"/>
			  <xs:element name="maxQueueSize" type="xs:integer"/>
			  <xs:element name="validationStringency" type="xs:string"/>
			</xs:sequence>
		  </xs:complexType>
		</xs:element>
		<xs:element name="samOptions">
		  <xs:complexType>
			<xs:sequence>
			  <xs:element name="picardBin" type="directoryPath"/>
			  <xs:element name="mergeLogBase" type="xs:integer"/>
			  <xs:element name="javaArgs" type="xs:string"/>
			  <xs:element name="qsubArgs" type="xs:string"/>
			  <xs:element name="validationStringency" type="xs:string"/>
			  <xs:element name="assumeSorted" type="xs:string"/>
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

	ValidatePath($data->{'srmaOptions'},         'srmaBinary',                               OPTIONAL); 
	ValidatePath($data->{'srmaOptions'},         'srmaJar',                                  OPTIONAL); 
	ValidatePath($data->{'srmaOptions'},         'javaBin',                                  OPTIONAL); 
	ValidateFile($data->{'srmaOptions'},         'javaArgs',                                 OPTIONAL);
	ValidatePath($data->{'srmaOptions'},         'qsubBin',                                  OPTIONAL); 
	ValidateOption($data->{'srmaOptions'},        'queueType',                                REQUIRED);

	# fix up queueType to allow the user NOT to specify 'arrayJob'
	# This is necessary because we're using XML::Simple for parsing
	if (ref(\$data->{'srmaOptions'}->{'queueType'}) eq "SCALAR") {
		my %hash = ();
		$hash{'content'} = $data->{'srmaOptions'}->{'queueType'};
		$hash{'arrayJob'} = "false";
		$data->{'srmaOptions'}->{'queueType'} = \%hash;
	}
	
	ValidateOption($data->{'srmaOptions'}->{'queueType'}, 'content',        		 REQUIRED);
	ValidateOptions($data->{'srmaOptions'}->{'queueType'}, 'content',       \%QUEUETYPES,	 REQUIRED);
	ValidateOption($data->{'srmaOptions'}->{'queueType'}, 'arrayJob',                        REQUIRED);

	ValidateOption($data->{'srmaOptions'},	     'offset',                                   OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'minimumAlleleProbability',		 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'minimumAlleleCoverage',			 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'maximumTotalCoverage',			 OPTIONAL);

	ValidateOption($data->{'srmaOptions'},       'referenceFasta',                           REQUIRED);
	ValidateOption($data->{'srmaOptions'}->{'referenceFasta'}, 'content',                    REQUIRED);
	ValidateFile($data->{'srmaOptions'}->{'referenceFasta'}, 'content',                      REQUIRED);
	ValidateOption($data->{'srmaOptions'}->{'referenceFasta'}, 'splitSize',                  REQUIRED);

	ValidatePath($data->{'srmaOptions'},         'runDirectory',                             REQUIRED); 
	ValidatePath($data->{'srmaOptions'},         'tmpDirectory',                             REQUIRED); 
	ValidateFile($data->{'srmaOptions'},         'inputBAMFile',				 REQUIRED);
	ValidateFile($data->{'srmaOptions'},         'outputBAMFile',				 REQUIRED);
	ValidateOption($data->{'srmaOptions'},         'range',			 		 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},         'ranges',				 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},         'validationStringency',			 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},         'correctBases',				 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},         'numThreads',				 OPTIONAL);
	ValidateOption($data->{'srmaOptions'},	     'maximumQueueSize',			 OPTIONAL);

	die "Attribute splitSize must be greater than or equal tozero.\n" if ($data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'} < 0);
	die "Attribute \"splitSize\" may not be used with \"ranges\".\n" if (defined($data->{'srmaOptions'}->{'ranges'}) && 0 < $data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'});
	die "One of 'srmaBinary' or 'srmaJar' must be present" if (!defined($data->{'srmaOptions'}->{'srmaJar'}) && !defined($data->{'srmaOptions'}->{'srmaBinary'}));
	die "Attribute 'srmaBinary' cannot be used when the 'srmaJar' attribute is present" if (defined($data->{'srmaOptions'}->{'srmaJar'}) && defined($data->{'srmaOptions'}->{'srmaBinary'}));
	die "Array jobs are only available for SGE queues.\n" if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'} && "true" eq $data->{'srmaOptions'}->{'queueType'}->{'arrayJob'});

	if(defined($data->{'srmaOptions'}->{'range'}) && defined($data->{'srmaOptions'}->{'ranges'})) {
		die("Both attributes \"range\" and \"ranges\" may not be used in conjuction.\n");
	}

	# picard
	if(defined($data->{'samOptions'})) {
		ValidatePath($data->{'samOptions'},       'picardBin',                                OPTIONAL); 
		ValidateOption($data->{'samOptions'},     'mergeLogBase',                             OPTIONAL);
		ValidateOption($data->{'samOptions'},     'cleanUpTmpDirectory',                      OPTIONAL);
		ValidateOption($data->{'samOptions'},     'qsubArgs',                                 OPTIONAL);
		ValidateOption($data->{'samOptions'},     'validationStringency',		      OPTIONAL);
		ValidateOption($data->{'samOptions'},     'assumeSorted',			      OPTIONAL);
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

	open(FH, "$fasta.dict") || die("Could not open FASTA dictionary: $fasta.dict");
	while(defined(my $line = <FH>)) {
		if($line =~ m/\@SQ\tSN:(\S+)\tLN:(\d+)\t/) {
			my @a = ();
			push(@a, $1); push(@a, $2);
			push(@$genomeInfo, \@a);
		}
	}
	close(FH);

	# TODO
	return 0;
}

sub CreateJobs {
	my ($data, $quiet, $start_step, $dryrun) = @_;

	my @srmaJobIDs = ();
	my @srmaOutputBAMs = ();
	my @samJobIDs = ();
	my @samOutputIDs = ();
	my %bamArrayJobs = ();

	# Create directories - Error checking...
	mkpath([$data->{'srmaOptions'}->{'runDirectory'}],    ($quiet) ? 0 : 1, 0755);
	mkpath([$data->{'srmaOptions'}->{'tmpDirectory'}],    ($quiet) ? 0 : 1, 0755);

	CreateJobsSRMA($data, $quiet, $start_step, $dryrun, \@srmaJobIDs, \@srmaOutputBAMs, \%bamArrayJobs);
	if(0 < scalar(@srmaJobIDs)) {
		CreateJobsSAM($data, $quiet, $start_step, $dryrun, \@srmaJobIDs, \@srmaOutputBAMs, \%bamArrayJobs);
	}

	if("PBS" ne $data->{'srmaOptions'}->{'queueType'}->{'content'}) {
		# remove holds
		foreach my $qsubID (@qsubGlobalIDs) {
			my $outID="";
			$outID=`qalter $qsubID -h U`; # SGE 
			chomp($outID);
			if($outID=~ m/Your job (\d+)/) {
				$outID= $1;
			}
			die unless (0 < length($outID));
			print STDERR "[srma submit] hold removed QSUBID=$qsubID\n";
		}
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
	my ($data, $quiet, $start_step, $dryrun, $qsubIDs, $outputBAMs, $bamArrayJobs) = @_;
	my @read_files = ();

	if($data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'} <= 0) { # do not split
		my $runFile = CreateRunFile($data, 'srma', "");
		my $cmd = "";

		if(defined($data->{'srmaOptions'}->{'srmaJar'})) { 
			$cmd = $data->{'srmaOptions'}->{'javaBin'} if defined($data->{'srmaOptions'}->{'javaBin'});
			$cmd .= "java";
			if(defined($data->{'srmaOptions'}->{'javaArgs'})) {
				$cmd .= " ".$data->{'srmaOptions'}->{'javaArgs'};
				if($data->{'srmaOptions'}->{'javaArgs'} !~ m/\-Xmx/) {
					$cmd .= " -Xmx2g";
				}
			}
			else {
				$cmd .= " -Xmx2g";
			}
			$cmd .= " -jar ".$data->{'srmaOptions'}->{'srmaJar'};
			$cmd .= " I=".$data->{'srmaOptions'}->{'inputBAMFile'};
			$cmd .= " O=".$data->{'srmaOptions'}->{'outputBAMFile'};
			$cmd .= " R=".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
			$cmd .= " RANGE=".$data->{'srmaOptions'}->{'range'} if(defined($data->{'srmaOptions'}->{'range'}));
			$cmd .= " RANGES=".$data->{'srmaOptions'}->{'ranges'} if(defined($data->{'srmaOptions'}->{'ranges'}));
			$cmd .= " OFFSET=".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
			$cmd .= " MINIMUM_ALLELE_PROBABILITY=".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
			$cmd .= " MINIMUM_ALLELE_COVERAGE=".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
			$cmd .= " MAXIMUM_TOTAL_COVERAGE=".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
			$cmd .= " CORRECT_BASES=".$data->{'srmaOptions'}->{'correctBases'} if(defined($data->{'srmaOptions'}->{'correctBases'}));
			$cmd .= " NUM_THREADS=".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
			$cmd .= " MAX_QUEUE_SIZE=".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
			$cmd .= " VALIDATION_STRINGENCY=".$data->{'srmaOptions'}->{'validationStringency'} if(defined($data->{'srmaOptions'}->{'validationStringency'}));
		}
		else {
			$cmd = $data->{'srmaOptions'}->{'srmaBinary'} if defined($data->{'srmaOptions'}->{'srmaBinary'});
			$cmd .= " -i ".$data->{'srmaOptions'}->{'inputBAMFile'};
			$cmd .= " -o ".$data->{'srmaOptions'}->{'outputBAMFile'};
			$cmd .= " -r ".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
			$cmd .= " -R ".$data->{'srmaOptions'}->{'range'} if(defined($data->{'srmaOptions'}->{'range'}));
			$cmd .= " -Z ".$data->{'srmaOptions'}->{'ranges'} if(defined($data->{'srmaOptions'}->{'ranges'}));
			$cmd .= " -O ".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
			$cmd .= " -p ".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
			$cmd .= " -c ".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
			$cmd .= " -t ".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
			$cmd .= " -C " if(defined($data->{'srmaOptions'}->{'correctBases'}));
			$cmd .= " -n ".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
			$cmd .= " -Q ".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
		}

		# Submit the job
		my @a = (); # empty array for job dependencies
		my $qsubID = SubmitJob($runFile, $quiet, ($start_step <= $STARTSTEP{"srma"}) ? 1 : 0, 0, $dryrun, $cmd, $data, 'srmaOptions', 'srma', \@a);
		# Do not add any SAM jobs (so no qsubIDs or outputBAMs)
	}
	else {
		my $start = 1;
		my $splitSize = $data->{'srmaOptions'}->{'referenceFasta'}->{'splitSize'};
		my @genomeInfo = ();

		getGenomeInfo($data->{'srmaOptions'}->{'referenceFasta'}->{'content'}, \@genomeInfo);

		my ($range_chr, $range_start, $range_end) = ("", -1, -1);
		if(defined($data->{'srmaOptions'}->{'range'})) {
			if($data->{'srmaOptions'}->{'range'} =~ m/(.+):(\d+)-(\d+)/) {
				$range_chr = $1;
				$range_start = $2;
				$range_end = $3;
				if($range_end < $range_start) {
					die("Range start position was greater than range end position");
				}
				my $found_chr = 0;
				for(my $i=0;$i<scalar(@genomeInfo);$i++) {
					if($range_chr eq $genomeInfo[$i]->[0]) {
						$found_chr = 1;
						last;
					}
				}
				if(0 == $found_chr) {
					die("Could not find chromosome [$range_chr] in the reference.\n");
				}
			}
			else {
				$range_chr = $data->{'srmaOptions'}->{'range'};
				my $found_chr = 0;
				for(my $i=0;$i<scalar(@genomeInfo);$i++) {
					if($range_chr eq $genomeInfo[$i]->[0]) {
						$range_start = 1;
						$range_end = $genomeInfo[$i]->[1];
						$found_chr = 1;
						last;
					}
				}
				if(0 == $found_chr) {
					die("Could not find chromosome [$range_chr] in the reference.\n");
				}
			}
		}

		my $finished_range = -1;
		for(my $i=0;$i<scalar(@genomeInfo);$i++) {
			my $chrName = $genomeInfo[$i]->[0];
			my $chrSize = $genomeInfo[$i]->[1];
			my $start = 1;
			if(defined($data->{'srmaOptions'}->{'range'})) {
				if($chrName ne $range_chr) {
					next;
				}
				$start = $range_start;
				$chrSize = $range_end;
			}

			# For SGE, submit an array job
			
			if ("true" eq $data->{'srmaOptions'}->{'queueType'}->{'arrayJob'}) {
				my $outputID = "$chrName";
				my $runFile = CreateRunFile($data, 'srma', $outputID);

				my $pre_cmd = "chrSize=$chrSize;\n";
				$pre_cmd .= "splitSize=$splitSize;\n";
				$pre_cmd .= "start=\$(( (\$SGE_TASK_ID - 1) * \$splitSize + 1));\n";
				$pre_cmd .= "end=\$((\$start+\$splitSize-1));\n";
				$pre_cmd .= "end=\$((\$chrSize < \$end ? \$chrSize : \$end));\n";
				$pre_cmd .= "outputFile=".$data->{'srmaOptions'}->{'tmpDirectory'}."srma.".$outputID."_\${start}-\${end}.bam;\n";

				my $cmd = "";
				if(defined($data->{'srmaOptions'}->{'srmaJar'})) { 
					$cmd .= $data->{'srmaOptions'}->{'javaBin'} if defined($data->{'srmaOptions'}->{'javaBin'});
					$cmd .= "java";
					if(defined($data->{'srmaOptions'}->{'javaArgs'})) {
						$cmd .= " ".$data->{'srmaOptions'}->{'javaArgs'};
						if($data->{'srmaOptions'}->{'javaArgs'} !~ m/-Xmx/) {
							$cmd .= " -Xmx2g";
						}
					}
					else {
						$cmd .= " -Xmx2g";
					}
					$cmd .= " -jar ".$data->{'srmaOptions'}->{'srmaJar'};
					$cmd .= " I=".$data->{'srmaOptions'}->{'inputBAMFile'};
					$cmd .= " O=\$outputFile";
					$cmd .= " R=".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
					$cmd .= " OFFSET=".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
					$cmd .= " MINIMUM_ALLELE_PROBABILITY=".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
					$cmd .= " MINIMUM_ALLELE_COVERAGE=".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
					$cmd .= " MAXIMUM_TOTAL_COVERAGE=".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
					$cmd .= " RANGE=\\\"$chrName\:\$start-\$end\\\"";
					$cmd .= " CORRECT_BASES=".$data->{'srmaOptions'}->{'correctBases'} if(defined($data->{'srmaOptions'}->{'correctBases'}));
					$cmd .= " NUM_THREADS=".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
					$cmd .= " MAX_QUEUE_SIZE=".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
					$cmd .= " VALIDATION_STRINGENCY=".$data->{'srmaOptions'}->{'validationStringency'} if(defined($data->{'srmaOptions'}->{'validationStringency'}));
				}
				else {
					$cmd = $data->{'srmaOptions'}->{'srmaBinary'} if defined($data->{'srmaOptions'}->{'srmaBinary'});
					$cmd .= " -i ".$data->{'srmaOptions'}->{'inputBAMFile'};
					$cmd .= " -o \$outputFile";
					$cmd .= " -r ".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
					$cmd .= " -O ".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
					$cmd .= " -p ".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
					$cmd .= " -c ".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
					$cmd .= " -t ".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
					$cmd .= " -R \\\"$chrName\:\$start-\$end\\\"";
					$cmd .= " -C " if(defined($data->{'srmaOptions'}->{'correctBases'}));
					$cmd .= " -n ".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
					$cmd .= " -q ".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
				}
				# Submit the job
				my @a = (); # empty array for job dependencies
				my $qsubID = SubmitArrayJob($runFile, $quiet, ($start_step <= $STARTSTEP{"srma"}) ? 1 : 0, 0, $dryrun, $pre_cmd, $cmd, $data, 
								 $start, ceil($chrSize/$splitSize), 1, 'srmaOptions', $outputID, \@a);
				push(@$qsubIDs, $qsubID) if (QSUBNOJOB ne $qsubID);

				for(;$start <= $chrSize;$start+=$splitSize) {
					my $end = $start + $splitSize - 1;
					if($chrSize < $end) { $end = $chrSize; }
					my $outputID = "$chrName\_$start\-$end";
					my $outputFile = CreateTmpOutputFile($data, 'srma', $outputID);
					push(@$outputBAMs, $outputFile);
					$bamArrayJobs->{$outputFile} = $qsubID;
				}
			}
			else {
				for(;$start <= $chrSize;$start+=$splitSize) {
					my $end = $start + $splitSize - 1;
					if($chrSize < $end) { $end = $chrSize; }
					my $within_range = -1;
					die("$start $end\n") unless ($start <= $end);
					my $outputID = "$chrName\_$start\-$end";
					my $runFile = CreateRunFile($data, 'srma', $outputID);
					my $outputFile = CreateTmpOutputFile($data, 'srma', $outputID);

					my $cmd = "";
					if(defined($data->{'srmaOptions'}->{'srmaJar'})) { 
						$cmd = $data->{'srmaOptions'}->{'javaBin'} if defined($data->{'srmaOptions'}->{'javaBin'});
						$cmd .= "java";
						if(defined($data->{'srmaOptions'}->{'javaArgs'})) {
							$cmd .= " ".$data->{'srmaOptions'}->{'javaArgs'};
							if($data->{'srmaOptions'}->{'javaArgs'} !~ m/-Xmx/) {
								$cmd .= " -Xmx2g";
							}
						}
						else {
							$cmd .= " -Xmx2g";
						}
						$cmd .= " -jar ".$data->{'srmaOptions'}->{'srmaJar'};
						$cmd .= " I=".$data->{'srmaOptions'}->{'inputBAMFile'};
						$cmd .= " O=$outputFile";
						$cmd .= " R=".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
						$cmd .= " OFFSET=".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
						$cmd .= " MINIMUM_ALLELE_PROBABILITY=".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
						$cmd .= " MINIMUM_ALLELE_COVERAGE=".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
						$cmd .= " MAXIMUM_TOTAL_COVERAGE=".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
						$cmd .= " RANGE=\\\"$chrName\:$start-$end\\\"";
						$cmd .= " CORRECT_BASES=".$data->{'srmaOptions'}->{'correctBases'} if(defined($data->{'srmaOptions'}->{'correctBases'}));
						$cmd .= " NUM_THREADS=".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
						$cmd .= " MAX_QUEUE_SIZE=".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
						$cmd .= " VALIDATION_STRINGENCY=".$data->{'srmaOptions'}->{'validationStringency'} if(defined($data->{'srmaOptions'}->{'validationStringency'}));
					}
					else {
						$cmd = $data->{'srmaOptions'}->{'srmaBinary'} if defined($data->{'srmaOptions'}->{'srmaBinary'});
						$cmd .= " -i ".$data->{'srmaOptions'}->{'inputBAMFile'};
						$cmd .= " -o $outputFile";
						$cmd .= " -r ".$data->{'srmaOptions'}->{'referenceFasta'}->{'content'};
						$cmd .= " -O ".$data->{'srmaOptions'}->{'offset'} if(defined($data->{'srmaOptions'}->{'offset'}));
						$cmd .= " -p ".$data->{'srmaOptions'}->{'minimumAlleleProbability'} if(defined($data->{'srmaOptions'}->{'minimumAlleleProbability'}));
						$cmd .= " -c ".$data->{'srmaOptions'}->{'minimumAlleleCoverage'} if(defined($data->{'srmaOptions'}->{'minimumAlleleCoverage'}));
						$cmd .= " -t ".$data->{'srmaOptions'}->{'maximumTotalCoverage'} if(defined($data->{'srmaOptions'}->{'maximumTotalCoverage'}));
						$cmd .= " -R \\\"$chrName\:$start-$end\\\"";
						$cmd .= " -C " if(defined($data->{'srmaOptions'}->{'correctBases'}));
						$cmd .= " -n ".$data->{'srmaOptions'}->{'numThreads'} if(defined($data->{'srmaOptions'}->{'numThreads'}));
						$cmd .= " -q ".$data->{'srmaOptions'}->{'maximumQueueSize'} if(defined($data->{'srmaOptions'}->{'maximumQueueSize'}));
					}

					# Submit the job
					my @a = (); # empty array for job dependencies
					my $qsubID = SubmitJob($runFile, $quiet, ($start_step <= $STARTSTEP{"srma"}) ? 1 : 0, 0, $dryrun, $cmd, $data, 'srmaOptions', $outputID, \@a);
					push(@$qsubIDs, $qsubID) if (QSUBNOJOB ne $qsubID);
					push(@$outputBAMs, $outputFile);
				}
			}
			if(0 == @$qsubIDs) {
				die;
			}
		}
	}
}


sub CreateJobsSAM {
	my ($data, $quiet, $start_step, $dryrun, $dependentQsubIDs, $dependentOutputBAMs, $bamArrayJobs) = @_;
	my @qsubIDs = ();
	my ($run_file, $outputID, $qsub_id);

	my $type = 'picard';
	my @outputBAMs = ();

	if(0 == scalar(@$dependentQsubIDs)) {
		# skip
		return;
	}

	if(!defined($data->{'samOptions'})) {
		# Delete the dependent jobs
		for(my $i=0;$i<scalar($dependentQsubIDs);$i++) {
			my $qsubID = $dependentQsubIDs->[$i];
			`qdel $qsubID`;
		}
		die("samOptions must be defined");
	}

	# Get dependent qsubIDs
	for(my $i=0;$i<scalar(@$dependentQsubIDs);$i++) {
		push(@qsubIDs, $dependentQsubIDs->[$i]);
	}

	# Get BAM file names
	for(my $i=0;$i<scalar(@$dependentOutputBAMs);$i++) {
		push(@outputBAMs, $dependentOutputBAMs->[$i]);
	}

	my $firstIter = 1;

	my $mergeLogBase = MERGE_LOG_BASE;
	if(defined($data->{'samOptions'}->{'mergeLogBase'}) && 1 < $data->{'samOptions'}->{'mergeLogBase'}) {
		$mergeLogBase = $data->{'samOptions'}->{'mergeLogBase'};
	}

	# Merge script(s)
	# Note: there could be too many dependencies, as well as too many output BAMs to merge.
	# Therefore, a hierarchical merge is required.
	my $mergeLevel = 0;
	my $shouldDepend = ($start_step <= $STARTSTEP{"srma"}) ? 1 : 0;
	do {
		$mergeLevel++;
		my $ctr = 0;
		my @curIDs = @qsubIDs;
		my @curBAMs = @outputBAMs;
		@qsubIDs = ();
		@outputBAMs = ();
		my $finalIter = 1;
		for(my $i=0;$i<scalar(@curBAMs);$i+=$mergeLogBase) {
			$ctr++;
		}
		if(1 < $ctr) {
			$finalIter = 0;
		}
		$ctr=0;
		for(my $i=0;$i<scalar(@curBAMs);$i+=$mergeLogBase) {
			$ctr++;
			# Get the subset of dependent jobs
			my @dependentIDs = ();
			my @dependentBAMs = ();
			for(my $j=$i;$j<scalar(@curBAMs) && $j<$i+$mergeLogBase;$j++) {
				push(@dependentBAMs, $curBAMs[$j]);
				push(@dependentIDs, $bamArrayJobs->{$curBAMs[$j]}) if ($firstIter && 
										       "true" eq $data->{'srmaOptions'}->{'queueType'}->{'arrayJob'});
				push(@dependentIDs, $curIDs[$j]) if (!$firstIter || "false" eq $data->{'srmaOptions'}->{'queueType'}->{'arrayJob'});
			}

			if ($firstIter && "true" eq $data->{'srmaOptions'}->{'queueType'}->{'arrayJob'}) {
				# The first time through, there will be many duplicate IDs if we used array jobs.
				# Remove them.
				my %depIDset = map { $_ => 1 } @dependentIDs;
				@dependentIDs = ();
				foreach my $id (%depIDset) {
					push (@dependentIDs, $id);
				}
			}

			# Set up ouptut ID, run file, and output BAM
			$outputID = "merge.$mergeLevel.$ctr";
			$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$type.".$outputID.".sh";
			my $outputBAM = "";
			if(0 == $finalIter) {
				$outputBAM = CreateTmpOutputFile($data, $type, $outputID);
			}
			else {
				$outputBAM = "".$data->{'srmaOptions'}->{'outputBAMFile'};
			}
			# Create the command
			if(!defined($data->{'samOptions'}->{'picardBin'})) { die("Picard bin required") };
			my $cmd = "";
			$cmd = $data->{'srmaOptions'}->{'javaBin'} if defined($data->{'srmaOptions'}->{'javaBin'});
			$cmd .= "java";
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
			foreach my $bam (@dependentBAMs) {
				$cmd .= " I=$bam";
			}
			$cmd .= " O=$outputBAM";
			$cmd .= " SO=coordinate";
			if(1 == $mergeLevel && defined($data->{'samOptions'}->{'assumeSorted'})) {
				# only on the first level, since it will then be sorted
				$cmd .= " AS=".$data->{'samOptions'}->{'assumeSorted'};
			}
			else {
				$cmd .= " AS=true";
			}
			$cmd .= " TMP_DIR=".$data->{'srmaOptions'}->{'tmpDirectory'};
			$cmd .= " VALIDATION_STRINGENCY=".$data->{'samOptions'}->{'validationStringency'} if(defined($data->{'samOptions'}->{'validationStringency'}));
			# Submit
			$qsub_id = SubmitJob($run_file, $quiet, ($start_step <= $STARTSTEP{"sam"}) ? 1 : 0, $shouldDepend, $dryrun, $cmd, $data, 'samOptions', $outputID, \@dependentIDs);
			if(QSUBNOJOB ne $qsub_id) {
				push(@qsubIDs, $qsub_id);
				push(@outputBAMs, $outputBAM);
			}
			else {
				# currently it must be submitted
				die;
			}
		}
		$shouldDepend = 1; # always depend on the next loops
		$firstIter = 0; # We are no longer on the first iteration
	} while(1 < scalar(@qsubIDs)); # while merging is necessary


# Clean up
	if(defined($data->{'srmaOptions'}->{'cleanUpTmpDirectory'}) && 1 == defined($data->{'srmaOptions'}->{'cleanUpTmpDirectory'})) {
		my @a = (); push(@a, $qsubIDs[0]); # push the merge/copy
		$outputID = "cleanup.tmpdirectory";
		$run_file = $data->{'srmaOptions'}->{'runDirectory'}."$outputID.sh";
		my $cmd = "rm -rv ".$data->{'srmaOptions'}->{'tmpDirectory'};
		$qsub_id = SubmitJob($run_file , $quiet, ($start_step <= $STARTSTEP{"sam"}) ? 1 : 0, 1, $dryrun, $cmd, $data, 'srmaOptions', $outputID, \@a);
	}
}

sub SubmitJob {
	my ($run_file, $quiet, $should_run, $should_depend, $dryrun, $command, $data, $type, $outputID, $dependent_jobIDs) = @_;
	if(0 < length($outputID)) {
		$outputID = "$type.$outputID"; $outputID =~ s/Options//g;
	}
	else {
		$outputID = "$type"; $outputID =~ s/Options//g;
	}

	if(!$quiet) {
		print STDERR "[srma submit] RUNFILE=$run_file\n";
	}

	if(1 == $should_run) {
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
		if(defined($data->{'srmaOptions'}->{'javaUse'})) {
			$output .= "use ".$data->{'srmaOptions'}->{'javaUse'}."\n";
		}
		$output .= "\nrun \"hostname\";\n";
		# Redirect PBS stderr/stdout, since it buffers them
		if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'}) {
			my $pbs_stderr_redirect = "$run_file.stderr.redirect";
			my $pbs_stdout_redirect = "$run_file.stdout.redirect";
			$output .= "run \"$command 2> $pbs_stderr_redirect > $pbs_stdout_redirect\";\n";
		}
		else {
			$output .= "run \"$command\";\n";
		}   
		$output .= "exit 0;\n";
		open(FH, ">$run_file") or die("Error.  Could not open $run_file for writing!\n");
		print FH "$output";
		close(FH);
	}

	# Create qsub command
	my $qsub = "";
	$qsub .= $data->{'srmaOptions'}->{'qsubBin'} if defined($data->{'srmaOptions'}->{'qsubBin'});
	$qsub .= "qsub" if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'}); # without a user hold
	$qsub .= "qsub -h" if ("SGE" eq $data->{'srmaOptions'}->{'queueType'}->{'content'}); # with a user hold, remove later

	if(0 < scalar(@$dependent_jobIDs) && 1 == $should_depend) {
		$qsub .= " -hold_jid ".join(",", @$dependent_jobIDs)         if ("SGE" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
		$qsub .= " -W depend=afterok:".join(":", @$dependent_jobIDs) if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
	}
	if(defined($data->{$type}->{'numThreads'}) && 1 < $data->{$type}->{'numThreads'}) {
		$qsub .= " -pe serial ".$data->{$type}->{'numThreads'}     if ("SGE" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
		$qsub .= " -l nodes=1:ppn=".$data->{$type}->{'numThreads'} if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
	}
	$qsub .= " ".$data->{$type}->{'qsubArgs'} if defined($data->{$type}->{'qsubArgs'});
	$qsub .= " -terse" if ("SGE" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
	$qsub .= " -N $outputID -o $run_file.out -e $run_file.err $run_file";

	if(1 == $should_run) {
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

		# save qsub id
		push(@qsubGlobalIDs, $qsub_id);

		if(!$quiet) {
			print STDERR "[srma submit] NAME=$outputID QSUBID=$qsub_id\n";
		}

		return $qsub_id;
	}
}

sub SubmitArrayJob {
	my ($run_file, $quiet, $should_run, $should_depend, $dryrun, $pre_command, $command, $data, 
	    $start, $end, $step, $type, $outputID, $dependent_jobIDs) = @_;

	if ("PBS" eq $data->{'srmaOptions'}->{'queueType'}->{'content'}) {
		die("Array jobs are only supported under SGE.");
	}

	if(0 < length($outputID)) {
		$outputID = "$type.$outputID"; $outputID =~ s/Options//g;
	}
	else {
		$outputID = "$type"; $outputID =~ s/Options//g;
	}

	if(!$quiet) {
		print STDERR "[srma submit] RUNFILE=$run_file\n";
	}

	if(1 == $should_run) {
		my $output = <<END_OUTPUT;
#\$ -e $run_file.\$TASK_ID.err
#\$ -o $run_file.\$TASK_ID.out
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
		if(defined($data->{'srmaOptions'}->{'javaUse'})) {
			$output .= "use ".$data->{'srmaOptions'}->{'javaUse'}."\n";
		}
		$output .= $pre_command;
		$output .= "\nrun \"hostname\";\n";
		$output .= "run \"$command\";\n";
		$output .= "exit 0;\n";
		open(FH, ">$run_file") or die("Error.  Could not open $run_file for writing!\n");
		print FH "$output";
		close(FH);
	}

	# Create qsub command
	my $qsub = "";
	$qsub .= $data->{'srmaOptions'}->{'qsubBin'} if defined($data->{'srmaOptions'}->{'qsubBin'});
	$qsub .= "qsub -h"; # with a user hold, remove later

	if(0 < scalar(@$dependent_jobIDs) && 1 == $should_depend) {
		$qsub .= " -hold_jid ".join(",", @$dependent_jobIDs);
	}
	if(defined($data->{$type}->{'numThreads'}) && 1 < $data->{$type}->{'numThreads'}) {
		$qsub .= " -pe serial ".$data->{$type}->{'numThreads'};
	}
	$qsub .= " ".$data->{$type}->{'qsubArgs'} if defined($data->{$type}->{'qsubArgs'});
	$qsub .= " -t ".$start."-".$end.":".$step;
	$qsub .= " -terse" if ("SGE" eq $data->{'srmaOptions'}->{'queueType'}->{'content'});
	$qsub .= " -N $outputID $run_file";

	if(1 == $should_run) {
		if(1 == $dryrun) {
			$FAKEQSUBID++;
			print STDERR "[srma submit] NAME=$outputID QSUBID=$FAKEQSUBID\n";
			return $FAKEQSUBID;
		}

		# Submit the qsub command
		my $qsub_id=`$qsub`;
		$qsub_id = "$qsub_id";
		chomp($qsub_id);

		if($qsub_id =~ m/(\d+)\..+/) {
			$qsub_id = $1;
		}
		die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
		if($qsub_id !~ m/^\S+$/) {
			die("Error submitting QSUB_COMMAND=$qsub\nQSUB_ID=$qsub_id\n") unless (0 < length($qsub_id));
		}

		# save qsub id
		push(@qsubGlobalIDs, $qsub_id);

		if(!$quiet) {
			print STDERR "[srma submit] NAME=$outputID QSUBID=$qsub_id\n";
		}

		return $qsub_id;
	}
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

=item B<-startstep>
Specifies on which step of the alignment process to start (default: srma). The
values can be "srma" or "sam".

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
