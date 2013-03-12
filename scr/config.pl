#!/usr/bin/perl -w

$flatlibpath = shift;
$flatlibpath =~ s#/$##;   # take away final / if any

my $n = "\n";

$file = "flatCommon.f90";
$fullfile = "$flatlibpath/src/".$file;


$template = &readFile($fullfile.".template");
my $progname = "subs";
my $ref_progname = \&$progname;
$myfile = $ref_progname->($template);
&printToFile($fullfile,$myfile); 


#readFile(path_to_file) -> turns back a string with the file
sub readFile { 
	my $path = shift;
	my $file = '';
	open (FILE,"< $path") or die "not able to open $path\.$n";
	while (<FILE>) {
		$file .= $_;
	}
	return $file;	
}


#printToFile(path_to_file,string)
sub printToFile {
	my $path = shift;
	my $file = shift;
	open (FILE,"> $path") or die "not able to open $path\.$n";
	print FILE $file;
}




sub subs {

   my $file = shift;

   my $start = "!Start path";
   my $end = "!End path";
   $str=" character (len=strlen), parameter :: flatlib_rootdir ='".$flatlibpath."/'".$n;
   $file =~ s/$start.+$end/$str/ms;

   return $file;

}



