#!/bin/bash

set +u
shopt -s extglob

######################## START BASIC FUNCTIONS ########################
# Number of arguments
nfield() { echo $# ; }

# Print the $1-th field of $2-${$#}
field() { local i=$1 ; shift ; ((i>0)) && echo "${@:$i:1}" ; }

swap_fields() {
 # Swap fields #$1 and #$2 in $3-$n
 local i1=$1 i2=$2 ix i1m i1p i2m i2p i12 i21
 shift 2
 if ((i1>0)) && ((i2>0)) && ((i1<=$#)) && ((i2<=$#)) && ((i1!=i2)) ; then
  if ((i2<i1)) ; then
   ix=$i1 ; i1=$i2 ; i2=$ix
  fi
  i1m=$((i1-1)) ; i1p=$((i1+1)) ; i2p=$((i2+1))
  i12=$((i2-i1-1))
  echo "${@:1:$i1m}" "${@:$i2:1}" "${@:$i1p:$i12}" "${@:$i1:1}" "${@:$i2p}"
 else
  echo "$@"
 fi
}

which_field() {
 # Output the field number of string $1 in $2-$n, or 0 if not present
 local i=0 str="$1"
 while (($#>1)) ; do i=$((i+1)) ; shift
  [ "$str" = "$1" ] && { echo $i ; return ; }
 done
 echo 0
}

in_line() {
 # Like which_field, but only set return value, no other output
 local str="$1"
 while (($#>1)) ; do shift ; [ "$str" = "$1" ] && return 0 ; done
 return 1
}

unpad() {
 # Remove leading and trailing blanks from "$@"
 local string="$@"
 while [ "${string:0:1}" = " " ] ; do
  string="${string:1}"
 done
 while [ "${string:$((${#string}-1)):1}" = " " ] ; do
  string="${string:0:$((${#string}-1))}"
 done
 echo "$string"
}

pretty_print() {
 # Print $3 with line folding at column $lwidth, with indentation $1 on the
 # first line and indentation $2 on the following.
 local indent1=$1 indent_rest=$2 text word line
 text="$(unpad "$3")"
 line=""
 while ((${#text}>0)) ; do
  word="${text%% *}"
  text="$(unpad "${text:${#word}}")"
  if [ -z "$line" ] ; then # only happens first time around
   line="$(printf "%${indent1}s")$word"
  else
   if ((${#line}+1+${#word}>lwidth)) ; then
    echo "$line$el" ; line="$(printf "%${indent_rest}s")$word"
   else
    line="$line $word"
   fi
  fi
 done
 [ -z "$line" ] || echo "$line$el"
 return 0
}

errstop() { echo "$el" ; pretty_print 1 3 "ERROR: $*" ; echo ; exit 1 ; }
errwarn() { echo -n "$el" ; pretty_print 1 3 "WARNING: $*" ; }

repeatl() {
 # Print $2 for the length of string $1
 local i=0 n=${#1} line=""
 while ((i<n)) ; do i=$((i+1)) ; line="$line$2" ; done
 echo "$line"
}

repeat() {
 # Print $2 $1 times
 local i=0 string=""
 while ((i<$1)) ; do i=$((i+1)) ; string="$string$2" ; done
 echo "$string"
}

# Add zeroes to the left of $2 until it has $1 characters.
zero_pad() { echo "$(repeat $(($1-${#2})) 0)$2" ; }

heading() {
 # Print $1 doubly underlined
 echo " $1$el"
 echo " $(repeatl "$1" =)$el"
}

lesser_heading() {
 # Print $1 underlined
 echo " $1$el"
 echo " $(repeatl "$1" -)$el"
}

# Print $1 spaces
space() { (($1>0)) && printf "%${1}s\n" "" ; }

check_number_N() { [[ "$1" == +([0-9]) ]] ; }

rel_versions() {
 # Print relationship (>/</=) between versions $1 and $2.
 local v1="$1" v2="$2" i1 i2
 while ((${#v1}>0)) && ((${#v2}>0)) ; do
  i1="${v1%%.*}" ; v1="${v1:${#i1}}" ; v1="${v1#.}"
  i2="${v2%%.*}" ; v2="${v2:${#i2}}" ; v2="${v2#.}"
  if check_number_N "$i1" && check_number_N "$i2" ; then
   ((i1>i2)) && { echo ">" ; return ; }
   ((i1<i2)) && { echo "<" ; return ; }
  else
   [[ "$i1" > "$i2" ]] && { echo ">" ; return ; }
   [[ "$i1" < "$i2" ]] && { echo "<" ; return ; }
  fi
 done
 ((${#v1}>0)) && { echo ">" ; return ; }
 ((${#v2}>0)) && { echo "<" ; return ; }
 echo "="
}

compare_versions() {
 # Compare version $1 against comma-separated list of patterns/version ranges
 # $2.
 local hvalue="$1" fvalue="$2" v
 [ -z "$hvalue" ] || [ -z "$fvalue" ] && return 20
 while ((${#fvalue}>0)) ; do
  v="${fvalue%%,*}"
  fvalue="$(unpad "${fvalue:${#v}}")"
  fvalue="$(unpad "${fvalue#,}")"
  v="$(unpad "$v")"
  compare_version "$hvalue" "$v" && return 0
 done
 return 1
}

compare_version() {
 # Compare version $1 against pattern/version range $2
 # Recognized forms of $2 are:
 # - '<bash-pattern>' (1 field); $1 is compared against pattern
 # - '<sign> <version>' (2 fields), where <sign> can be <, <=, >=, >; $1
 #   is compared against <version> to see if <sign> correctly describes their
 #   relationship
 # - '<version1> - <version2>' (3 fields); $1 is compared against <version1>
 #   and <version2> to check that <version1> <= $1 <= <version2>
 local version="$1" pattern="$2" rel1 rel2
 set -- $pattern
 case $# in
 1) [[ "$version" == $pattern ]] && return 0 ;;
 2)
  rel1=$(rel_versions $version $2)
  case "$1:$rel1" in
  ">=:>"|">=:="|"<=:<"|"<=:="|">:>"|"<:<") return 0 ;;
  esac ;;
 3)
  rel1=$(rel_versions $2 version)
  rel2=$(rel_versions $version $3)
  case "$rel1:$rel2" in
  "<:<"|"=:<"|"<:="|"=:=") return 0 ;;
  esac ;;
 esac
 return 1
}

progress_bar() {
 # Print a progress bar for action $1 at point $2 of $3.  No terminal codes
 # on output, that should be handled by the calling routine.
 local barwidth=50 title=$1 i=$2 n=$3 pc pc1 pc2 isp b1 b2 s
 local spinner="|/-\\"
 pc=" $(((100*i)/n))% "
 pc2=$((barwidth/2+3))
 pc1=$((pc2-${#pc}+1))
 isp=$(((barwidth*i)/n))
 s="${spinner:$((i%${#spinner})):1}"
 if ((isp==0)) ; then
  bar="$s$(space $((pc1-1)))$pc$(space $((barwidth-pc2)))| $title"
 elif ((isp<pc1)) ; then
  bar="|$(repeat $((isp-1)) =)$s$(space $((pc1-isp-1)))$pc$(space\
   $((barwidth-pc2)))| $title"
 elif ((isp<=pc2)) ; then
  bar="|$(repeat $((pc1-1)) =)$pc$(space $((barwidth-pc2)))| $title"
 elif ((isp<barwidth)) ; then
  bar="|$(repeat $((pc1-1)) =)$pc$(repeat $((isp-pc2)) =)$s$(space\
   $((barwidth-isp-1)))| $title"
 else
  bar="|$(repeat $((pc1-1)) =)$pc$(repeat $((barwidth-pc2)) =)| $title"
 fi
 echo "$bar"
}

mktemp_aix() {
mkdir -p .casinotemp_$$
echo .casinotemp_$$
}

######################### END BASIC FUNCTIONS #########################

########################### START FUNCTIONS ###########################
usage_and_exit() {
 # Display error message $1 (optionally), then display usage and exit
 cat <<_EOH
${0##*/}: helper script to pick a suitable CASINO_ARCH

Usage:
 ./${0##*/} [<options>]

Options for operation mode selection:
 --help
   Display this help and exit.

 --list  --list-ext  --list-all
   List generic/extended/all available CASINO_ARCHs.

 --auto
   List all CASINO_ARCHs matching the current machine.

 --gen
   Generate a new CASINO_ARCH for the current machine - this is an
   interactive mode.

Options for all modes:
  --type=<type>
    Restrict the inspected CASINO_ARCHs to those of TYPE=<type>, where <type>
    can be one of "single", "parallel" and "cluster".

  --dump
    Do not use terminal codes for "fancy" output; instead write plain text
    suitable for redirecting the output to a file.

Options in --list* mode:
  --details
    Show details for each CASINO_ARCH

Options in --auto mode:
  --explain
    Print the reasons why CASINO_ARCHs are accepted as matches for the current
    machine.

  --explain-all
    Print the reasons why CASINO_ARCHs are accepted or rejected as matches for
    the current machine.

  --nfork=<nfork>
    Spawn <nfork> forks to parse the .arch files in --auto mode.  By default
    <nfork> is twice the number of cores of the host machine.
_EOH
 (($#>0)) && errstop "$1"
 exit 0
}

initialize() {
 local maxproc proc_per_fork
 # Set script name.
 script_name=${0##*/}
 # Select temporary directory.  This is just for placing named pipes.
 tmpdir="$TMPDIR" ; [ -z "$tmpdir" ] && tmpdir="$HOME"
 # Define tags
 define_tags
 # Define a set of tags to issue warnings about
 mandatory_taglist="TYPE F90 CC"
 recommended_taglist="DESCRIPTION MAINTAINER DATE ARCH KERNEL F90_VERSION\
  COMMAND_CHECK_F90_VERSION COMMAND_CHECK_F90"
 # Define substitutable variables in this script
 variable_list="F90 CC CXX"
 # Defaults
 explain=0 ; fancy=1 ; span=all ; details=0
 debug_compiler_output=/dev/null
 # Set number of forks to the number of (virtual) cores.
 nfork=1
 [ -e "/proc/cpuinfo" ] && nfork=$(cat /proc/cpuinfo | grep -cE "^[pP]rocessor")
 check_number_N "$nfork" || nfork=1
 # Reduce number of forks according to limit on number of user processes.
 proc_per_fork=10 # over-estimate of number of processes required to run a fork
 maxproc=$(ulimit -u 2> /dev/null)
 if check_number_N "$maxproc" ; then
  ((nfork*proc_per_fork>maxproc)) && nfork=$((maxproc/proc_per_fork))
 fi
 ((nfork<1)) && nfork=1
 # mktemp command to create a temporary directory (on $TMPDIR or /tmp).
 case "$(mktemp --version 2>&1)" in
 *GNU*) # GNU mktemp
  mktemp_d_command="mktemp -d" ;;
 *found*) # AIX, or possibly others where mktemp doesn't exist
  mktemp_d_command="mktemp_aix" ;;
 *) # assume BSD's mktemp (as found in *BSD and Mac OS)
  mktemp_d_command="mktemp -d -t tmp" ;;
 esac
}

parse_cmd_line() {
 # Read command line.
 run_mode=list ; target_TYPE="" ; span=generic
 while (($#>0)) ; do
  case "$1" in
  --help) usage_and_exit ;;
  --list) run_mode=list ; span=generic ;;
  --list-all) run_mode=list ; span=all ;;
  --list-ext) run_mode=list ; span=ext ;;
  --auto) run_mode=auto ; span=all ;;
  --gen) run_mode=gen ;;
  --explain) explain=1 ;;
  --explain-all) explain=2 ;;
  --details) details=1 ;;
  --type=*) target_TYPE=${1#--type=} ;;
  --dump) fancy=0 ;;
  --debug-compiler) debug_compiler_output=/dev/stdout ;;
  --nfork=*)
   nfork=${1#--nfork=}
   check_number_N $nfork && ((nfork>0)) || usage_and_exit "In \
    --nfork=<nfork>, <nfork> must be a positive integer." ;;
  *) usage_and_exit "Unknown option '$1'." ;;
  esac
  shift
 done
 [ "$run_mode" = list ] || span=all
}

post_config() {
 if ((fancy==0)) ; then
  el="" ; cr=""
 elif [ -z "$el" ] || [ -z "$cr" ] ; then
  el="" ; cr="" ; fancy=0
 fi
}

build_file_list() {
 # Build list of files
 local CASINO_ARCH CASINO_ARCH_GENERIC
 nfile=0 ; filelist=""
 for file in *.arch ; do
  [ -d "$file" ] && continue
  CASINO_ARCH=${file%.arch}
  CASINO_ARCH_GENERIC=${CASINO_ARCH%.*}
  case "$span" in
  generic) [ "$CASINO_ARCH" = "$CASINO_ARCH_GENERIC" ] || continue ;;
  ext) [ "$CASINO_ARCH" = "$CASINO_ARCH_GENERIC" ] && continue ;;
  esac
  filelist="$filelist $file"
  nfile=$((nfile+1))
 done
}

init_tput() {
 # Get terminal codes
 local IFS_save="$IFS" term_columns
 [[ "$TERM" == xterm-* ]] && export TERM=xterm
 if tput -S < /dev/null >& /dev/null ; then
  # Quick version
  {
   IFS=$(echo -e "\t")
   read el cr term_columns
   IFS="$IFS_save"
  } < <(echo -e "el \nht\n cr \nht\n cols \n" | tput -S)
  (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 else
  el=$(tput el 2>/dev/null)
  cr=$(tput cr 2>/dev/null)
  term_columns=$(tput cols)
 fi
 [ -z "$el" ] && el=$(tput ce 2>/dev/null)
 check_number_N "$term_columns" || term_columns=80
 lwidth=$((term_columns-1))
}

get_host_params() {
 # Detect current machine parameters.
 host_ARCH=$(uname -m 2> /dev/null)
 host_KERNEL=$(uname -s 2> /dev/null)
 host_DISTRIBUTION=$(
  if type -P lsb_release >& /dev/null ; then
   lsb_release -ds 2> /dev/null
  elif [ -e /etc/redhat-release ] ; then
   head -n 1 /etc/redhat-release 2> /dev/null
  elif [ -e /etc/SuSE-release ] ; then
   head -n 1 /etc/SuSE-release 2> /dev/null
  elif [ -e /etc/slackware-version ] ; then
   head -n 1 /etc/slackware-version 2> /dev/null
  elif [ -e /etc/gentoo-release ] ; then
   head -n 1 /etc/gentoo-release 2> /dev/null
  elif [ -e /etc/arch-release ] ; then
   echo "Arch Linux"
  elif [ -e /etc/debian_version ] ; then
   echo "Debian GNU/Linux $(head -n 1 /etc/debian_version 2> /dev/null)"
  fi
 )
 host_DISTRIBUTION=${host_DISTRIBUTION#\"}
 host_DISTRIBUTION=${host_DISTRIBUTION%\"}
 host_OS=$(uname -o 2> /dev/null)
 host_HOSTNAME=$(hostname 2> /dev/null)
 host_DOMAIN=$(hostname -d 2> /dev/null)
}

match_host_param() {
 # Match parameter $1 against value for this machine.
 local var="$1" hvalue fvalue v found
 eval hvalue=\"\$host_$var\"
 eval fvalue=\"\$$var\"
 [ -z "$fvalue" ] || [ -z "$hvalue" ] && return 20
 found=0
 while ((${#fvalue}>0)) ; do
  v="${fvalue%%,*}"
  fvalue="$(unpad "${fvalue:${#v}}")"
  fvalue="$(unpad "${fvalue#,}")"
  v="$(unpad "$v")"
  [[ "$hvalue" == $v ]] && return 0
 done
 return 1
}

autodetect_arch() {
 local ifork pipe flist nfile nfile_per_fork ffiles c s r nfork_actual
# local list_match_host_distro list_match_host list_match_distro
# local list_match_other
 local t1 t2 nfork_extra nfile_per_fork_extra n any_printed
 local no_output=0
 list_match_host_distro=""
 list_match_host=""
 list_match_distro=""
 list_match_other=""
 while (($#>0)) ; do
  case $1 in
  --no-output) no_output=1 ;;
  esac
  shift
 done
 # Print heading when required.
 case "$explain" in
 0) : ;;
 1) heading "Detailed match criteria for accepted CASINO_ARCHs" ;;
 2) heading "Detailed match criteria for all CASINO_ARCHs" ;;
 esac
 # Unsort the filelist to balance the load among forks.
 flist=$(\
  {
   for f in $filelist ; do
    echo $f
   done
  } | { \
   while read line ; do \
    echo "$(zero_pad 5 $RANDOM)$line" ; \
   done ; \
  } | sort | cut -b6- \
 )
 nfile=$(nfield $flist)
 # Set up fork parameters.
 nfork_actual=$nfork
 ((nfile<nfork_actual)) && nfork_actual=$nfile
 nfile_per_fork=$((nfile/nfork_actual))
 nfork_extra=$((nfile%nfork_actual))
 nfile_per_fork_extra=$((nfile_per_fork+1))
 # Create named pipe
 pipe="$tmpdir/.tmp_${script_name}"
 [ -e "$pipe" ] && rm -f "$pipe"
 [[ "$(uname -s)" != CYGWIN* ]] && mkfifo "$pipe" # Cygwin won't do this
 # Loop over forks, handing each of them a part of the file list.
 ifork=0 ; while ((ifork<nfork_actual)) ; do ifork=$((ifork+1))
  if ((fancy==1)) ; then
   t1=" $(progress_bar "Preparing" $ifork $nfork_actual)"
   echo -n "$t1$el$cr"
  fi
  if ((ifork<=nfork_extra)) ; then
   n=$nfile_per_fork_extra
  else
   n=$nfile_per_fork
  fi
  ffiles="$(set -- $flist ; echo "${@:1:$n}")"
  flist="$(set -- $flist ; echo "${@:$((n+1))}")"
  if [[ "$(uname -s)" != CYGWIN* ]] ; then # asynchronous processing
   sleep 1 >> "$pipe" &
   autodetect_arch_fork "$ffiles" >> "$pipe" &
  else # synchronous processing (Cygwin only)
   autodetect_arch_fork "$ffiles" >> "$pipe" 
  fi
 done
 ((fancy==1)) && echo -n " $(progress_bar "Searching" 0 $nfile) $el$cr"
 # Receive data from forks.
 {
  ifile=0 ; while read c s r ; do ifile=$((ifile+1))
   if ((fancy==1)) ; then
    t1=$(nfield $list_match_host_distro $list_match_host $list_match_distro\
     $list_match_other)
    if ((t1==0)) ; then
     t2=" $(progress_bar "Searching" $ifile $nfile)"
    elif ((t1==1)) ; then
     t2=" $(progress_bar "Searching ($t1 match)" $ifile $nfile)"
    else
     t2=" $(progress_bar "Searching ($t1 matches)" $ifile $nfile)"
    fi
    echo -n "$t2$el$cr"
   fi
   case $s in
   0)
    if ((explain>1)) ; then
     case "$(set -- $r ; echo $1)" in
     F90|CC|CXX)
      pretty_print 1 3 "X $c NOT compatible due to missing executable\
       '$(set -- $r ; echo $2)'" ;;
     *_CHECK)
      pretty_print 1 3 "X $c NOT compatible due to correctly-named compiler\
      '$(set -- $r ; echo $2)' printing a version string which does not\
      match known patterns for the required compiler" ;;
     *_VERSION)
       pretty_print 1 3 "X $c NOT compatible due to compiler version\
        '$(set -- $r ; echo $2)' not matching '$(set -- $r ; echo ${@:3})'" ;;
     *) pretty_print 1 3 "X $c NOT compatible due to $r mismatch" ;;
     esac
     ((fancy==1)) && echo -n "$t2$el$cr"
    fi ;;
   1)
    if in_line HOSTNAME $r || in_line DOMAIN $r ; then
     if in_line DISTRIBUTION $r ; then
      list_match_host_distro="$list_match_host_distro $c"
     else
      list_match_host="$list_match_host $c"
     fi
    else
     if in_line DISTRIBUTION $r ; then
      list_match_distro="$list_match_distro $c"
     else
      list_match_other="$list_match_other $c"
     fi
    fi
    if ((explain>0)) && [ ! -z "$r" ] ; then
     pretty_print 1 3 "* $c matches all of $r"
     ((fancy==1)) && echo -n "$t2$el$cr"
    fi ;;
   esac
  done
 } < $pipe
 echo -n "$el"
 rm -f "$pipe"
 ((ifile!=nfile)) && errwarn "Received $ifile reports instead of $nfile?"
 # Sort matches alphabetically and print them.
 # (NOTE: in sed command use literal newlines not \n or e.g. BSD/AIX 
 # won't work.. Also require single quotes not double.
 [ -z "$list_match_host_distro" ] || list_match_host_distro=$(echo\
  $list_match_host_distro | sed 's/ /\
/g' | sort)
 [ -z "$list_match_host" ] || list_match_host=$(echo $list_match_host\
  | sed 's/ /\
/g' | sort)
 [ -z "$list_match_distro" ] || list_match_distro=$(echo $list_match_distro\
  | sed 's/ /\
/g' | sort)
 [ -z "$list_match_other" ] || list_match_other=$(echo $list_match_other\
  | sed 's/ /\
/g' | sort)
 ((explain>0)) && echo
 if ((no_output==0)) ; then
  any_printed=0
  if [ ! -z "$list_match_host_distro" ] ; then
   ((any_printed==1)) && echo
   heading "CASINO_ARCHs matching this machine's hostname and Linux\
  distribution"
   for CASINO_ARCH in $list_match_host_distro ; do
    pretty_print 1 3 "* $CASINO_ARCH"
   done
   any_printed=1
  fi
  if [ ! -z "$list_match_host" ] ; then
   ((any_printed==1)) && echo
   heading "CASINO_ARCHs matching this machine's hostname"
   for CASINO_ARCH in $list_match_host ; do
    pretty_print 1 3 "* $CASINO_ARCH"
   done
   any_printed=1
  fi
  if [ ! -z "$list_match_distro" ] ; then
   ((any_printed==1)) && echo
   heading "CASINO_ARCHs matching this machine's Linux distribution"
   for CASINO_ARCH in $list_match_distro ; do
    pretty_print 1 3 "* $CASINO_ARCH"
   done
   any_printed=1
  fi
  if [ ! -z "$list_match_other" ] ; then
   ((any_printed==1)) && echo
   heading "CASINO_ARCHs matching this machine's basic parameters"
    for CASINO_ARCH in $list_match_other ; do
    pretty_print 1 3 "* $CASINO_ARCH"
   done
   any_printed=1
  fi
 fi
}

autodetect_arch_fork() {
 # Match current architecture against CASINO_ARCHs in list $1, and produce
 # a single-line report for each of them, for parsing in the main
 # autodetection function above.
 local flist="$1" file match_var tlist_s tlist_b
 for file in $flist ; do
  CASINO_ARCH=${file%.arch}
  CASINO_ARCH_GENERIC=${CASINO_ARCH%.*}
  match_var=""
  clear_tags
  # Test for TYPE first.
  if [ ! -z "$target_TYPE" ] ; then
   quickload_tags_make "TYPE" "$file"
   if [ "$target_TYPE" != "$TYPE" ] ; then
    echo "$CASINO_ARCH 0 TYPE"
    continue
   fi
  fi
  # Test basic tags.
  quickload_tags_scalar "HOSTNAME DOMAIN" "$file"
  set -- $(check_host)
  (($1==0)) && echo "$CASINO_ARCH 0 ${@:2}" && continue
  match_var="$match_var ${@:2}"
  quickload_tags_scalar "ARCH KERNEL DISTRIBUTION OS" "$file"
  set -- $(check_OS)
  (($1==0)) && echo "$CASINO_ARCH 0 ${@:2}" && continue
  match_var="$match_var ${@:2}"
  # Test for existence of compilers.
  quickload_tags_make "F90 CC CXX ENVIRONMENT_COMMAND" "$file"
  [ -z "$ENVIRONMENT_COMMAND" ] && ENVIRONMENT_COMMAND=:
  var_F90="$ENVIRONMENT_COMMAND ; $F90"
  var_CC="$ENVIRONMENT_COMMAND ; $CC"
  var_CXX="$ENVIRONMENT_COMMAND ; $CXX"
  set -- $(check_compilers F90 CC CXX)
  (($1==0)) && echo "$CASINO_ARCH 0 ${@:2}" && continue
  tlist_s="${@:3:$2}" ; tlist_b="${@:$(($2+3))}"
  # Deep-test compilers.
  if [ ! -z "$tlist_s" ] ; then
   quickload_tags_scalar "$tlist_s" "$file"
   quickload_tags_block "$tlist_b" "$file"
   set -- $(check_compilers_deep F90 CC CXX)
   (($1==0)) && echo "$CASINO_ARCH 0 ${@:2}" && continue
   match_var="$match_var ${@:2}"
  fi
  # Return outcome
  echo "$CASINO_ARCH 1 $match_var"
 done
}

check_host() {
 # Test basic tags.
 # Must have loaded: HOSTNAME
 local var match_var=""
 if [ -z "$DOMAIN" ] ; then
  for var in HOSTNAME ; do
   match_host_param $var
   case $? in
   0) match_var="$match_var $var" ;;
   1) echo "0 $var" ; return 1 ;;
   esac
  done
 else # some distinct machines have identical HOSTNAMEs but different DOMAINs
  for var in HOSTNAME DOMAIN ; do
   match_host_param $var
   case $? in
   0) match_var="$match_var $var" ;;
   1) echo "0 $var" ; return 1 ;;
   esac
  done
 fi
 echo "1 $match_var"
 return 0
}

check_OS() {
 # Test basic tags.
 # Must have loaded: ARCH KERNEL DISTRIBUTION OS
 local var match_var=""
 for var in ARCH KERNEL DISTRIBUTION OS ; do
  match_host_param $var
  case $? in
  0) match_var="$match_var $var" ;;
  1) echo "0 $var" ; return 1 ;;
  esac
 done
 echo "1 $match_var"
 return 0
}

check_compilers() {
 # Test for existance of compilers $*.
 # Must have loaded: F90 CC CXX ENVIRONMENT_COMMAND
 local var match_var="" tlist_s tlist_b t1 c fvalue
 tlist_s="" ; tlist_b=""
 for var in $* ; do
  eval "fvalue=\"\$$var\""
  [ -z "$fvalue" ] && continue
  c=$(\
   set -- $fvalue ;\
   while (($#>0)) ; do\
    [[ "$1" == *=* ]] || break ;\
    shift ;\
   done ;\
   echo $1\
  )
  t1=$(
   bash < <(
    echo "$ENVIRONMENT_COMMAND"
    echo "type -P \"$c\" >& /dev/null && echo 1 || echo 0"
   ) 2> /dev/null
   (exit 0) # bypass bash fd leak (v3.2 - v4.1)
  )
  if ((t1==0)) ; then
   echo "0 $var $c"
   return 1
  fi
  tlist_s="$tlist_s ${var}_VERSION"
  tlist_b="$tlist_b COMMAND_CHECK_$var COMMAND_CHECK_${var}_VERSION"
 done
 echo "1 $(nfield $tlist_s) $tlist_s $tlist_b"
}

check_compilers_deep() {
 # Deep-test compilers $*.
 # Must have loaded: output from check_compilers
 local var size fvalue version host_version match_var
 match_var=""
 # Match compilers
 for var in $* ; do
  # Match compiler identification
  eval "size=\"\$COMMAND_CHECK_$var\""
  if [ ! -z "$size" ] ; then
   eval_tag COMMAND_CHECK_$var
   if [ "$(eval_code_tmp COMMAND_CHECK_${var})" != 1 ] ; then
    eval "fvalue=\"\$$var\""
    echo "0 ${var}_CHECK $fvalue"
    return 1
   fi
   match_var="$match_var ${var}_CHECK"
  fi
  # Match compiler version
  eval "version=\"\$${var}_VERSION\" ;\
   size=\"\$COMMAND_CHECK_${var}_VERSION\""
  if [ ! -z "$version" ] && [ ! -z "$size" ] ; then
   eval_tag COMMAND_CHECK_${var}_VERSION
   host_version="$(eval_code_tmp COMMAND_CHECK_${var}_VERSION)"
   if [ -z "$host_version" ] ; then
    echo "0 ${var}_VERSION"
    return 1
   fi
   eval "version=\"\$${var}_VERSION\""
   compare_versions "$host_version" "$version"
   case "$?" in
   0) match_var="$match_var ${var}_VERSION" ;;
   1) echo "0 ${var}_VERSION $host_version $version"
    return 1 ;;
   esac
  fi
 done
 echo "1 $match_var"
}

list_arch() {
 # List available architectures in class $1 (all/generic/ext).
 local span="$1"
 local ifile file CASINO_ARCH CASINO_ARCH_GENERIC tag nlines ivector t1
 if [ ! -z "$target_TYPE" ] ; then
  case "$span" in
  all)
   heading "List of all supported CASINO_ARCHs (TYPE=$target_TYPE)" ;;
  generic)
   heading "List of supported generic CASINO_ARCHs (TYPE=$target_TYPE)" ;;
  ext)
   heading "List of supported extended CASINO_ARCHs (TYPE=$target_TYPE)" ;;
  esac
 else
  case "$span" in
  all) heading "List of all supported CASINO_ARCHs" ;;
  generic) heading "List of supported generic CASINO_ARCHs" ;;
  ext) heading "List of supported extended CASINO_ARCHs" ;;
  esac
 fi
 # Loop over files and print information.
 ifile=0 ; for file in $filelist ; do ifile=$((ifile+1))
  if ((fancy==1)) ; then
   t1=" $(progress_bar "Reading" $ifile $nfile)"
   echo -n "$t1$el$cr"
  fi
  CASINO_ARCH=${file%.arch}
  CASINO_ARCH_GENERIC=${CASINO_ARCH%.*}
  if ((details==0)) && [ -z "$target_TYPE" ] ; then
   pretty_print 1 3 "* $CASINO_ARCH"
   ((fancy==1)) && echo -n "$t1$el$cr"
  else
   clear_tags
   # Skip unwanted types
   if [ ! -z "$target_TYPE" ] ; then
    quickload_tags_make TYPE "$file"
    [ "$target_TYPE" != "$TYPE" ] && continue
   fi
   # Print
   if ((details==1)) ; then
    if [ -z "$target_TYPE" ] ; then
     quickload_tags_scalar "DESCRIPTION MAINTAINER DATE COMMENT TYPE\
      QUEUEING_SYSTEM ARCH KERNEL OS HOSTNAME DOMAIN" "$file"
    else
     quickload_tags_scalar "DESCRIPTION MAINTAINER DATE COMMENT\
      QUEUEING_SYSTEM ARCH KERNEL OS HOSTNAME DOMAIN" "$file"
    fi
    heading $CASINO_ARCH
    for tag in DESCRIPTION MAINTAINER DATE COMMENT TYPE QUEUEING_SYSTEM ARCH\
     KERNEL OS HOSTNAME DOMAIN ; do
     eval "[ -z \"\$$tag\" ] || pretty_print 1 3 \"$tag: \$$tag\""
     ((fancy==1)) && echo -n "$t1$el$cr"
    done
    echo "$el"
    ((fancy==1)) && echo -n "$t1$el$cr"
   else
    pretty_print 1 3 "* $CASINO_ARCH"
    ((fancy==1)) && echo -n "$t1$el$cr"
   fi
  fi
 done
 ((fancy==1)) && echo -n "$cr$el$cr"
}

prompt() {
 local prompt_arg prompt_var prompt_prompt prompt_default prompt_allow
 local prompt_val prompt_non_empty prompt_force_integer prompt_remember
 local prev_value what_default interpret_equal=0
 for prompt_arg in "$@" ; do
  case "$prompt_arg" in
  *=*) eval prompt_${prompt_arg%%=*}=\"\${prompt_arg#*=}\" ;;
  *) prompt_var="$prompt_arg" ;;
  esac
  shift
 done
 [ -z "$prompt_var" ] && return 1
 [ -z "$prompt_prompt" ] && prompt_prompt=CHOICE
 [ -z "$prompt_remember" ] && prompt_remember=0
 what_default="DEFAULT"
 if ((prompt_remember==1)) ; then
  eval "prev_value=\"\$$prompt_var\""
  if [ ! -z "$prev_value" ] ; then
   prompt_default="$prev_value"
   what_default="PREVIOUSLY-ENTERED"
  fi
 fi
 while : ; do
  if [ -z "$prompt_default" ] ; then
   read -e -p " $prompt_prompt> " $prompt_var || return 127
  else
   case "$BASH_VERSION" in
   0.*|1.*|2.*|3.*)
    echo " TYPE '=' FOR $what_default VALUE: $prompt_default"
    read -e -p " $prompt_prompt> " $prompt_var || return 127
    interpret_equal=1  ;;
   4.*)
    read -e -p " $prompt_prompt> " -i "$prompt_default" $prompt_var \
     || return 127 ;;
   esac
  fi
  echo
  [ "$prompt_non_empty" = 1 ] || [ ! -z "$prompt_allow" ]\
   || [ ! -z "$prompt_force_integer" ] || ((interpret_equal==1))\
   && eval prompt_val=\"\$$prompt_var\"
  ((interpret_equal==1)) && [ "$prompt_val" = "=" ]\
   && eval "$prompt_var=\"\$prompt_default\""\
   && prompt_val="$prompt_default"
  [ "$prompt_non_empty" = 1 ] && [ -z "$prompt_val" ]\
   && echo " Must provide a non-empty value, try again." && continue
  [ "$prompt_force_integer" = 1 ] && ! check_number_N "$prompt_val"\
   && echo " Must provide an integer, try again." && continue
  [ ! -z "$prompt_allow" ] && ! in_line "$prompt_val" $prompt_allow\
   && echo " Wrong option, try again." && continue
  return 0
 done
}

gen_arch() {
 # Guide the user through the process of creating a new CASINO_ARCH
 local file line os_match cut_here_string cluster_variables
 local cluster_variables_internal cluster_variables_user var desc deflt
 local chosen_os chosen_f90 chosen_cc chosen_cxx chosen_compiler chosen_queue
 local chosen_f90_native chosen_cc_native
 local suggest_TYPE suggest_choice suggest_CPN suggest_arch
 local arch_form review_compiler review_queue arch_file cmt1 cmt
 local tlist_s_f90 tlist_b_f90 tlist_s_cc tlist_b_cc tlist_s_cxx tlist_b_cxx
 local tlist_s_f90_native tlist_b_f90_native
 local tlist_s_cc_native tlist_b_cc_native
 local tlist_s tlist_b ignore_bintest bin_needs_mpirun
 local ignore_mpitest ignore_mpitest_requested
 local f90_is_for f90_var_name f90_var_val cc_is_for cc_var_name cc_var_val
 local c2fflags_var_name c2fflags_var_val env_cmd
 local tmp_dir t1 choice i c cc clow cclow chuman allow bscript skip comp_desc
 local queue last_ncore os comp chosen pref_editor name email
 local -a match amatch match_ldpath vmatch
 local user_TYPE user_ENVIRONMENT_COMMAND user_F90 user_CC user_CXX
 local user_ENVIRONMENT_COMMAND_NATIVE
 local user_F90_NATIVE user_CC_NATIVE user_CFLAGS_F90_INTERFACE_NATIVE
 local user_FFLAGS_opt_NATIVE
 local user_NEED_ETIME user_INCLUDE_DIR user_LIB_PATH user_LDLIBS_all
 local user_CFLAGS_F90_INTERFACE user_CFLAGS_ETIME user_FFLAGS_opt
 local user_FFLAGS_debug user_FFLAGS0_libs
 local user_HAVE_BLAS user_HAVE_LAPACK
 local user_FFLAGS_OPENMP_yes user_CFLAGS_OPENMP_yes user_SUPPORT_SHM
 local user_CFLAGS_SHM
 local user_CORES_PER_NODE user_RUN_PARALLEL user_CORES_PER_NODE_CLUSTER
 local user_QUEUEING_SYSTEM user_SUBMIT_SCRIPT user_MAX_NCORE user_MAX_WALLTIME
 local user_MAX_NCORE_per_queue user_MAX_WALLTIME_per_queue
 local user_MAX_NCORE_interdep user_MAX_WALLTIME_interdep
 local user_HOSTNAME user_CASINO_ARCH user_DESCRIPTION user_MAINTAINER
 local user_UTILS_MODE
 local user_DATE user_MPI_VERSION
 local -a user_SCRIPT_HEAD user_SCRIPT_RUN
 local IFS_save="$IFS"
 local stages_done next_stage
 local stage_name_TYPE="Type of machine"
 local stage_name_XCOMP="Cluster cross-compilation"
 local stage_name_STAGING="Multiple file systems and file staging"
 local stage_name_ENV="Environment configuration"
 local stage_name_COMPILERS="Compiler configuration"
 local stage_name_TOPOLOGY="Runtime configuration"
 local stage_name_ARCH_DETAILS="Additional information"
 # Load 'prompt' locally so that it works on local variables
 eval $(type prompt | tail -n +2)
 # Internal functions
 stage_heading() {
  local stage_name
  eval "stage_name=\"\$stage_name_$current_stage\""
  [ -z "$stage_name" ] && stage_name="Unknown stage"
  heading "$stage_name"
 }
 step_heading() {
  lesser_heading "$1"
 }
 meta_menu() {
  local action stage i n allow_list stage_name
  echo ; echo
  pretty_print 1 1 "--Entering meta-menu--"
  echo
  pretty_print 1 1 "Choose an option:"
  allow_list=""
  i=0 ; for stage in $stages_done ; do i=$((i+1))
   eval "stage_name=\"\$stage_name_$stage\""
   pretty_print 1 5 "[$i] Go back to: $stage_name"
   allow_list="$allow_list $i"
  done
  n=$i
  eval "stage_name=\"\$stage_name_$current_stage\""
  pretty_print 1 5 "[r] Restart current stage: $stage_name"
  allow_list="$allow_list r"
  pretty_print 1 5 "[a] Abort guided configuration"
  allow_list="$allow_list a"
  pretty_print 1 5 "[c] Continue guided configuration"
  allow_list="$allow_list c"
  while ! prompt action allow="$allow_list" ; do echo ; done
  pretty_print 1 1 "--Exiting meta-menu--"
  echo
  case "$action" in
  a) return 1 ;;
  c) return 0 ;;
  r) next_stage=$current_stage  ; return 2 ;;
  *)
   if check_number_N "$action" ; then
    set -- $stages_done
    next_stage="${@:$action:1}"
    stages_done="${@:1:$((action-1))}"
    return 2
   fi ;;
  esac
  return 0
 }
 # Print header
 cat <<_EOF

 ================================
 GUIDED CASINO_ARCH CONFIGURATION
 ================================

_EOF
 # Print important information
 pretty_print 1 12 "IMPORTANT: you will need to provide information about\
  your machine, like compiler names, required environment settings, job\
  limits, submission scripts (where applicable), etc.  This set-up tool will\
  be as helpful as possible but it does not do magic - make sure you have the\
  machine documentation at hand."
 echo
 pretty_print 1 12 "IMPORTANT: at any point in the set-up process you may\
  want to go back to an earlier step or abort the process entirely.  You can\
  press Ctrl+D at any prompt (having deleted any existing text after the\
  prompt) to bring up a meta-menu that will offer these options."
 echo
 # Detect OS from existing definitions
 heading "OS detection [automatic]"
 echo -n " Detecting OS...$el$cr"
 for file in os/* ; do
  os=${file%.arch}
  os=${os##*/}
  clear_tags
  # Load OS tags.
  quickload_tags_scalar "ARCH KERNEL DISTRIBUTION OS" "$file"
  set -- $(check_OS)
  (($1==0)) && continue
  os_match="$os_match $os"
 done
 echo -n "$el$cr"
 case "$(nfield $os_match)" in
 1)
  chosen_os="$(unpad $os_match)"
  clear_tags
  quickload_tags_scalar "DESCRIPTION" "os/$chosen_os.arch"
  pretty_print 1 1 "Detected OS: $DESCRIPTION" ;;
 *)
  chosen_os=""
  pretty_print 1 1 "No predefined OS detected." ;;
 esac
 echo

 # Loop over stages.
 next_stage=TYPE
 while : ; do

  current_stage="$next_stage"
  next_stage=""
  [ "$current_stage" != WRITE ] && stage_heading

  case "$current_stage" in

  TYPE)
   # Ask user to provide type
   suggest_TYPE=""
   if type -P qstat >& /dev/null || type -P showq >& /dev/null\
    || type -P qsub >& /dev/null || type -P bsub >& /dev/null\
    || type -P llsubmit >& /dev/null || type -P llq >& /dev/null ; then
    suggest_TYPE=cluster
   elif [ ! -z "$chosen_os" ] ; then
    clear_tags
    load_tags "os/$chosen_os.arch"
    eval_tag CORES_PER_NODE
    check_number_N $CORES_PER_NODE || CORES_PER_NODE=1
    ((CORES_PER_NODE==1)) && suggest_TYPE=single || suggest_TYPE=parallel
   else
    suggest_TYPE=single
   fi
   pretty_print 1 1 "What type of machine is this?  At the prompt below enter:"
   pretty_print 1 3 "- 'single' if this is a single-processor (single-core)\
    workstation"
   pretty_print 1 3 "- 'parallel' if this is a multi-processor (multi-core)\
    workstation"
   pretty_print 1 3 "- 'cluster' if this is a cluster, that is, a collection\
    of single- or multi-processor nodes, each of which runs an independent\
    copy of the operating system, and which usually has a queueing system to\
    run computational jobs on them"
   while : ; do
    prompt user_TYPE prompt=TYPE default="$suggest_TYPE"\
     allow="single parallel cluster" remember=1 && break
    meta_menu
    case $? in
    0) : ;;
    1) return 1 ;;
    2) continue 2 ;;
    esac
   done

   if [ "$user_TYPE" = cluster ] ; then
    next_stage=XCOMP
   else
    next_stage=ENV
   fi ;;

  XCOMP)
   # Ask whether a cross-compiler is needed
   user_UTILS_MODE=0
   pretty_print 1 1 "The compute nodes of some clusters have CPUs that differ\
    from those in the login nodes, and different compilers are required\
    to produce binaries that work on each of them - namely the main CASINO\
    binary runs on the compute nodes and the CASINO utilities run on the\
    login nodes.  The compiler that produces binaries for the compute nodes\
    is called a cross compiler, while that whose binaries work on the login\
    nodes is a regular native compiler."
   echo
   pretty_print 1 1 "Most often clusters don't require cross-compilation; if\
    your cluster documentation does not mention any of this, the option to\
    choose below is likely to be number 2."
   echo
   pretty_print 1 1 "Choose an option:"
   pretty_print 1 5 "[1] This machine uses a cross-compiler"
   pretty_print 1 5 "[2] This machine does not use a cross-compiler"
   while : ; do
    prompt choice allow="1 2" default=2 && break
    meta_menu
    case $? in
    0) : ;;
    1) return 1 ;;
    2) continue 2 ;;
    esac
   done
   case "$choice" in
   1) user_UTILS_MODE=native ;;
   2) : ;;
   esac

   next_stage=STAGING ;;

   STAGING)
   # Ask about separate file systems and staging
   pretty_print 1 1 "The compute nodes of some clusters have completely\
   separate file systems from the login nodes (i.e. they are essentially\
   different computers). Such machines require\
   'staging' (physical transfer of the CASINO executable and input/output\
   files back and forth between the two file systems at run-time) and other\
   procedures which the CASINO run script needs to know about."
   echo
   pretty_print 1 1 "Such setups are very rare (Japan's K computer being\
   one example); for almost all other machines you should choose option 2\
   below."
   echo
   pretty_print 1 1 "Choose an option:"
   pretty_print 1 5 "[1] This machine requires file staging"
   pretty_print 1 5 "[2] This machine does not require file staging"
   while : ; do
    prompt choice allow="1 2" default=2 && break
    meta_menu
    case $? in
    0) : ;;
    1) return 1 ;;
    2) continue 2 ;;
    esac
   done
   case "$choice" in
   1) user_RELPATHNAMES="yes" ;;
   2) : ;;
   esac

   next_stage=ENV ;;

  ENV)
   # Ask the user to configure the environment
   pretty_print 1 1 "On some machines, you need to run a series of commands\
    before being able to use the Fortran compiler.  You should enter these\
    commands below."
   echo
   pretty_print 1 1 "Specifically, you need this if your machine requires you\
    to:"
   echo
   pretty_print 1 3 "* modify the PATH or LD_LIBRARY_PATH environment\
    variables"
   pretty_print 3 3 "Example:"
   pretty_print 5 7 "export PATH=/f90/bin:\${PATH}\
    LD_LIBRARY_PATH=/f90/lib:\${LD_LIBRARY_PATH}"
   echo
   pretty_print 1 3 "* use 'module load' commands"
   pretty_print 3 3 "Example:"
   pretty_print 5 7 ". /etc/profile.d/modules.sh ; module purge ; module load\
    default-impi"
   pretty_print 3 3 "NB, we recommend sourcing the modules.sh file (whatever\
    its name and location) explicitly, preferrably followed by a 'module\
    purge', and only then load the specific modules you need, so that changes\
    to the default module set do not break your CASINO set-up."
   echo
   pretty_print 1 3 "* set specific environment variables"
   pretty_print 3 3 "Example:"
   pretty_print 5 7 "export MPI=OpenMPI TCM_IFORT_VER=11.1-011"
   echo
   pretty_print 1 1 "If you do not need any of the above, leave the answer\
    below empty."
   echo
   pretty_print 1 1 "Notes:"
   pretty_print 1 3 "- Use curly braces to reference environment variables:\
    \${PATH}"
   pretty_print 1 3 "- Use back-quotes for command substitution: \`uname -i\`"
   pretty_print 1 3 "- Use semicolons ';' to separate multiple commands"
   echo
   if [ "$user_UTILS_MODE" != native ] ; then
    pretty_print 1 1 "Enter the required environment-modifying command, or\
     empty to omit:"
    while : ; do
     prompt user_ENVIRONMENT_COMMAND prompt=ENVIRONMENT_COMMAND remember=1\
      && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    if [ -z "$user_ENVIRONMENT_COMMAND" ] ; then
     user_ENVIRONMENT_COMMAND=:
     pretty_print 1 1 "ENVIRONMENT_COMMAND left unset."
     echo
    fi
   else
    step_heading "Cross-compiler environment"
    pretty_print 1 1 "Enter the required environment-modifying command for\
     cross-compilation and for running cross-compiled jobs on the cluster\
     compute nodes, or empty to omit:"
    while : ; do
     prompt user_ENVIRONMENT_COMMAND prompt=ENVIRONMENT_COMMAND remember=1\
      && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    if [ -z "$user_ENVIRONMENT_COMMAND" ] ; then
     user_ENVIRONMENT_COMMAND=:
     pretty_print 1 1 "ENVIRONMENT_COMMAND for cross-compiler left unset."
     echo
    fi
    step_heading "Native compiler environment"
    pretty_print 1 1 "Enter the required environment-modifying command for\
     native compilation and for running native binaries on the login nodes,\
     or empty to omit:"
    while : ; do
     prompt user_ENVIRONMENT_COMMAND_NATIVE prompt=ENVIRONMENT_COMMAND\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    if [ -z "$user_ENVIRONMENT_COMMAND_NATIVE" ] ; then
     user_ENVIRONMENT_COMMAND_NATIVE=:
     pretty_print 1 1 "ENVIRONMENT_COMMAND for native compiler left unset."
     echo
    fi
   fi

   next_stage=COMPILERS ;;

  COMPILERS)
   # Ask the user to configure the compilers
   if [ "$user_UTILS_MODE" = native ] ; then
    pretty_print 1 1 "We will now configure both the Fortran and (possibly)\
     C compilers.  For each of them you will need to provide the name of a\
     cross-compiler (the one used for producing binaries to be run on the\
     compute nodes of this cluster) and the native compiler (used for\
     light-weight utilities which run on the log-in node)."
    echo
   fi
   ignore_bintest_requested=0 ; bin_needs_mpirun=0
   for f90_is_for in src utils ; do
    eval "ignore_bintest=\"\$ignore_bintest_requested_$f90_is_for\""
    case "$f90_is_for" in
    src)
     f90_var_name=user_F90
     f90_var_val=""
     env_cmd="$user_ENVIRONMENT_COMMAND"
     if [ "$user_UTILS_MODE" != native ] ; then
      step_heading "Fortran compiler name"
     else
      step_heading "Fortran cross-compiler name"
      ignore_bintest=1
     fi ;;
    utils)
     [ "$user_UTILS_MODE" != native ] && continue
     f90_var_name=user_F90_NATIVE
     f90_var_val=""
     env_cmd="$user_ENVIRONMENT_COMMAND_NATIVE"
     step_heading "Native Fortran compiler name" ;;
    esac
    while : ; do
     if [ "$user_UTILS_MODE" = native ] ; then
      case "$f90_is_for" in
      src)
       pretty_print 1 1 "Enter the name of the Fortran 95 cross-compiler\
        binary - ideally a wrapper that automatically includes the MPI\
        libraries (e.g. 'mpif90'), possibly followed by any required\
        architecture selection options (e.g. '-target=barcelona'):" ;;
      utils)
       pretty_print 1 1 "Enter the name of the Fortran 95 native compiler\
        binary (e.g. 'gfortran'), possibly followed by any required\
        architecture selection options (e.g. '-target=native'):" ;;
      esac
     else
      if [ "$user_TYPE" = single ] ; then
       pretty_print 1 1 "Enter the name of the Fortran 95 compiler binary (e.g.\
        'gfortran'):"
      else
       pretty_print 1 1 "Enter the name of the Fortran 95 compiler binary -\
        ideally a wrapper that automatically includes the MPI libraries (e.g.\
        'mpif90'):"
      fi
     fi
     while : ; do
      prompt $f90_var_name prompt=F90 non_empty=1 remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 4 ;;
      esac
     done
     eval "f90_var_val=\"\$$f90_var_name\""
     case "$f90_is_for" in
     src) tlist_s_f90="" ; tlist_b_f90="" ;;
     utils) tlist_s_f90_native="" ; tlist_b_f90_native="" ;;
     esac
     # Test compiler
     pretty_print 1 1 "Testing Fortran compiler:"
     # Check compiler exists
     clear_tags
     F90="$f90_var_val"
     ENVIRONMENT_COMMAND="$env_cmd"
     var_F90="$ENVIRONMENT_COMMAND ; $f90_var_val"
     set -- $(check_compilers F90)
     if (($1==0)) ; then
      pretty_print 1 1 "Could not find compiler, try again. If you need to\
       load 'modules' or modify your PATH press Ctrl+D at an empty prompt\
       to go back to the '$stage_name_ENV' stage."
      echo ; continue
     fi
     case "$f90_is_for" in
     src)
      tlist_s_f90="${@:3:$2}" ; tlist_b_f90="${@:$(($2+3))}" ;;
     utils)
      tlist_s_f90_native="${@:3:$2}"
      tlist_b_f90_native="${@:$(($2+3))}" ;;
     esac
     pretty_print 1 1 "- Fortran compiler exists"
     # Test compilation
     run_compiler_test compile f90="$F90" env="$ENVIRONMENT_COMMAND"\
      output=1 mpirun="$user_RUN_PARALLEL" use_mpirun="$bin_needs_mpirun"
     t1=$?
     ((ignore_bintest==1)) && ((t1>1)) && t1=0
     case $t1 in
     0) : ;;
     1)
      pretty_print 1 1 "This does not appear to be a working Fortran 95\
       compiler.  If you need to load 'modules' or modify environment\
       variables press Ctrl+D at an empty prompt to go back to the\
       '$stage_name_ENV' stage."
      continue ;;
     2|3)
      pretty_print 1 1 "This Fortran compiler produces binaries which do not\
       appear to run properly."
      pretty_print 1 1 "Choose an option:"
      allow_list=""
      if [ "$user_TYPE" != single ] && ((bin_needs_mpirun==0)) ; then
       pretty_print 1 5 "[m] The MPI library requires the use of 'mpirun' (or\
        its equivalent) to run binaries linked against it (e.g. HP-MPI is\
        known to require this)"
       allow_list="$allow_list m"
      fi
      if [ "$user_TYPE" = cluster ] && [ "$user_UTILS_MODE" != native ] ; then
       pretty_print 1 5 "[x] This machine uses a cross-compiler; I would like\
        to go back to the '$stage_name_XCOMP' stage"
       allow_list="$allow_list x"
      fi
      if [ "$user_TYPE" = cluster ] && [ "$user_UTILS_MODE" != native ] ; then
       pretty_print 1 5 "[e] This compiler requires setting environment\
        variables to work properly; I would like to go back to the
        '$stage_name_ENV' stage"
       allow_list="$allow_list e"
      fi
      pretty_print 1 5 "[i] I would like to ignore this error"
      allow_list="$allow_list i"
      pretty_print 1 5 "[f] I would like to type the Fortran compiler name\
       again"
      allow_list="$allow_list f"
      while : ; do
       prompt choice allow="$allow_list" && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 4 ;;
       esac
      done
      case "$choice" in
      m)
       bin_needs_mpirun=1
       pretty_print 1 1 "Enter the command required to run an MPI calculation.\
        Here you must use the following variables:"
       pretty_print 1 3 "- &NPROC& expands to the number of MPI processes"
       pretty_print 1 3 "- &BINARY& expands to the name of the binary"
       while : ; do
        while : ; do
         prompt user_RUN_PARALLEL prompt=RUN_PARALLEL\
          default="mpirun -np &NPROC& &BINARY&" non_empty=1 remember=1 && break
         meta_menu
         case $? in
         0) : ;;
         1) return 1 ;;
         2) continue 5 ;;
         esac
        done
        ( [[ "$user_RUN_PARALLEL" != *"&NPROC&"* ]]\
         || [[ "$user_RUN_PARALLEL" != *"&BINARY&"* ]] )\
         && pretty_print 1 1 "The command must contain the &NPROC& and &BINARY&\
         variables, try again." && echo && continue
        break
       done
       continue ;;
      x)
       t1=""
       for stage in $stages_done ; do
        [ "$stage" = XCOMP ] && break
        t1="$t1 $stage"
       done
       stages_done="$t1" ; next_stage=XCOMP
       continue 3 ;;
      e)
       t1=""
       for stage in $stages_done ; do
        [ "$stage" = ENV ] && break
        t1="$t1 $stage"
       done
       stages_done="$t1" ; next_stage=ENV
       continue 3 ;;
      i) eval "ignore_bintest_requested_$f90_is_for=1" ; ignore_bintest=1 ;;
      f) continue ;;
      esac ;;
     esac
     ((ignore_bintest==0)) && pretty_print 1 1 "- Fortran compiler seems to\
      work."
     # Test for ETIME extension
     if [ "$f90_is_for" = src ] ; then
      run_compiler_test internal_etime f90="$F90"\
       env="$ENVIRONMENT_COMMAND" output=1 mpirun="$user_RUN_PARALLEL"\
       use_mpirun="$bin_needs_mpirun"
      t1=$?
      ((ignore_bintest==1)) && ((t1>1)) && t1=0
      case "$t1" in
      1)
       user_NEED_ETIME=yes
       pretty_print 1 3 "- Fortran compiler needs external ETIME: C compiler\
        required." ;;
      *)
       user_NEED_ETIME=no
       pretty_print 1 3 "- Fortran compiler has internal ETIME extension." ;;
      esac
     fi
     # Check for automatic inclusion of BLAS/LAPACK
     run_compiler_test blas_lapack f90="$F90"\
      env="$ENVIRONMENT_COMMAND" output=1 mpirun="$user_RUN_PARALLEL"\
      use_mpirun="$bin_needs_mpirun"
     t1=$?
     ((t1>1)) && ((ignore_bintest==1)) && t1=0
     if [ "$user_UTILS_MODE" != native ] ; then
      case "$t1" in
      0)
       pretty_print 1 3 "- Fortran compiler includes BLAS/LAPACK"
       user_HAVE_BLAS=yes
       user_HAVE_LAPACK=yes ;;
      *)
       pretty_print 1 3 "- Fortran compiler does not include BLAS/LAPACK"
       user_HAVE_BLAS=no
       user_HAVE_LAPACK=no ;;
      esac
     else
      case "$t1.$f90_is_for" in
      0.src)
       pretty_print 1 3 "- Fortran cross-compiler includes BLAS/LAPACK"
       user_HAVE_BLAS=yes
       user_HAVE_LAPACK=yes ;;
      0.utils)
       pretty_print 1 3 "- Native Fortran compiler includes BLAS/LAPACK"
       user_HAVE_BLAS_NATIVE=yes
       user_HAVE_LAPACK_NATIVE=yes ;;
      *.src)
       pretty_print 1 3 "- Fortran cross-compiler does not include BLAS/LAPACK"
       user_HAVE_BLAS=no
       user_HAVE_LAPACK=no ;;
      *.utils)
       pretty_print 1 3 "- Native Fortran compiler does not include\
        BLAS/LAPACK"
       user_HAVE_BLAS_NATIVE=no
       user_HAVE_LAPACK_NATIVE=no ;;
      esac
     fi
     echo
     # MPI libraries
     if [ "$user_TYPE" != single ] && [ "$f90_is_for" = src ] ; then
      step_heading "MPI library"
      user_MPI_VERSION=2
      user_INCLUDE_DIR="" ; user_LIB_PATH="" ; user_LDLIBS_all=""
      # Test MPI
      while : ; do
       run_compiler_test mpi f90="$F90"\
        fflags="$user_INCLUDE_DIR" lflags="$user_LIB_PATH $user_LDLIBS_all"\
        env="$ENVIRONMENT_COMMAND" output=2\
        mpirun="$user_RUN_PARALLEL" use_mpirun="$bin_needs_mpirun"
       t1=$?
       ((ignore_bintest==1)) && ((t1>1)) && { user_MPI_VERSION=1 ; t1=0 ; }
       ((t1==2)) && { user_MPI_VERSION=1 ; t1=0 ; }
       case "$t1" in
       1|3)
        echo
        pretty_print 1 3 "The Fortran compiler does not correctly include\
         and/or link an MPI library."
        echo
        pretty_print 1 1 "Choose an option:"
        pretty_print 1 5 "[l] I would like to enter the -I, -L and -l options\
         manually"
        pretty_print 1 5 "[f] I have entered the wrong Fortran compiler name;\
         I would like to restart the '$stage_name_COMPILERS' stage"
        pretty_print 1 5 "[i] I would like to ignore this error"
        while : ; do
         prompt choice allow="l f i" && break
         meta_menu
         case $? in
         0) : ;;
         1) return 1 ;;
         2) continue 4 ;;
         esac
        done
        case "$choice" in
        l)
         pretty_print 1 1 "Enter any -I options required to include the MPI\
          libraries, for example '-I/opt/OpenMPI/include' (the directory after\
          -I should be the location of the mpif.h header file):"
         while : ; do
          prompt user_INCLUDE_DIR prompt=INCLUDE_DIR remember=1 && break
          meta_menu
          case $? in
          0) : ;;
          1) return 1 ;;
          2) continue 4 ;;
          esac
         done
         pretty_print 1 1 "Enter any -L options required to link the MPI\
          libraries, for example '-L/opt/OpenMPI/lib' (the directory after -L\
          should be the location of the libmpi.a library file, or of any\
          otherwise named files providing the MPI library):"
         while : ; do
          prompt user_LIB_PATH prompt=LIB_PATH remember=1 && break
          meta_menu
          case $? in
          0) : ;;
          1) return 1 ;;
          2) continue 4 ;;
          esac
         done
         pretty_print 1 1 "Enter any -l options required to link the MPI\
          libraries, for example '-lmpi' (-l<name> means that the library file\
          in the above location is called 'lib<name>.a'):"
         while : ; do
          prompt user_LDLIBS_all prompt=LDLIBS remember=1 && break
          meta_menu
          case $? in
          0) : ;;
          1) return 1 ;;
          2) continue 4 ;;
          esac
         done ;;
        f)
         t1=""
         for stage in $stages_done ; do
          [ "$stage" = COMPILERS ] && break
          t1="$t1 $stage"
         done
         stages_done="$t1" ; next_stage=COMPILERS
         continue 3 ;;
        i) ignore_mpitest_requested=1 ; ignore_mpitest=1 ; break ;;
        esac ;;
       0)
        if [ "$user_UTILS_MODE" = native ] ; then
         pretty_print 1 3 "- Fortran compiler correctly includes and links an\
          MPI library (MPI version cannot be tested; will not use MPI 2\
          features to ensure compatibility)"
        else
         pretty_print 1 3 "- Fortran compiler correctly includes and links an\
          MPI $user_MPI_VERSION library."
        fi
        echo
        break ;;
       esac
      done
     fi
     break
    done
   done
   for cc_is_for in src utils ; do
    cc_compulsory=0
    case "$cc_is_for" in
    src)
     [ "$user_NEED_ETIME" = yes ] && cc_compulsory=1
     f90_var_name=user_F90
     f90_var_val="$user_F90"
     cc_var_name=user_CC
     c2fflags_var_name=user_CFLAGS_F90_INTERFACE
     c2fflags_var_val=""
     env_cmd="$user_ENVIRONMENT_COMMAND"
     if [ "$user_UTILS_MODE" != native ] ; then
      step_heading "C compiler name"
     else
      step_heading "C cross-compiler name"
     fi ;;
    utils)
     [ "$user_UTILS_MODE" != native ] && continue
     f90_var_name=user_F90_NATIVE
     f90_var_val="$user_F90_NATIVE"
     cc_var_name=user_CC_NATIVE
     cc_var_val=""
     c2fflags_var_name=user_CFLAGS_F90_INTERFACE_NATIVE
     c2fflags_var_val=""
     env_cmd="$user_ENVIRONMENT_COMMAND_NATIVE"
     step_heading "Native C compiler name" ;;
    esac
    while : ; do
     eval "ignore_bintest=\"\$ignore_bintest_requested_$cc_is_for\""
     if [ "$user_UTILS_MODE" = native ] ; then
      case "$cc_is_for" in
      src)
       if [ "$user_NEED_ETIME" = yes ] ; then
        pretty_print 1 1 "A C cross-compiler is required because the selected\
         Fortran compiler does not provide the ETIME extension.  Specifying a\
         C cross-compiler also allows you to build:"
        pretty_print 1 3 "- the optional \"shared-memory blips\" feature of\
         CASINO, for which an MPI wrapper around the C compiler should be\
         provided"
        echo
        pretty_print 1 1 "Enter the name of the C cross-compiler binary (e.g.\
         'mpicc'), possibly followed by any required architecture selection\
         options (e.g. '-target=barcelona'):"
       else
        pretty_print 1 1 "A C cross-compiler is optional.  Specifying a C\
         cross-compiler allows you to build:"
        pretty_print 1 3 "- the optional \"shared-memory blips\" feature of\
         CASINO, for which an MPI wrapper around the C compiler should be\
         provided"
        echo
        pretty_print 1 1 "Enter the name of the C cross-compiler binary (e.g.\
         'mpicc'), possibly followed by any required architecture selection\
         options (e.g. '-target=barcelona'), or empty to omit:"
       fi
       ignore_bintest=1 ;;
      utils)
       pretty_print 1 1 "A native C compiler is optional.  Specifying a native\
        C compiler allows you to build:"
       pretty_print 1 3 "- the optional \"jeep_to_pwfn\" JEEP-to-CASINO\
        wave-function converter"
       echo
       pretty_print 1 1 "Enter the name of the native C compiler binary (e.g.\
        'gcc'), possibly followed by any required architecture selection\
        options (e.g. '-target=native'):"
      esac
     else
      if [ "$user_NEED_ETIME" = yes ] ; then
       pretty_print 1 1 "A C compiler is required because the selected Fortran\
        compiler does not provide the ETIME extension.  Specifying a C\
        compiler also allows you to build:"
       if [ "$user_TYPE" != single ] ; then
        pretty_print 1 3 "- the optional \"shared-memory blips\" feature of\
         CASINO, for which an MPI wrapper around the C compiler should be\
         provided"
       fi
       pretty_print 1 3 "- the optional \"jeep_to_pwfn\" JEEP-to-CASINO\
        wave-function converter"
       echo
       pretty_print 1 1 "Enter the name of the C compiler binary (e.g. 'gcc'):"
      else
       pretty_print 1 1 "A C compiler is optional.  Specifying a C compiler\
        allows you to build:"
       if [ "$user_TYPE" != single ] ; then
        pretty_print 1 3 "- the optional \"shared-memory blips\" feature of\
         CASINO, for which an MPI wrapper around the C compiler should be\
         provided"
       fi
       pretty_print 1 3 "- the optional \"jeep_to_pwfn\" JEEP-to-CASINO\
        wave-function converter"
       echo
       pretty_print 1 1 "Enter the name of the C compiler binary (e.g. 'gcc'),\
        or empty to omit:"
      fi
     fi
     while : ; do
      prompt $cc_var_name prompt=CC non_empty=$cc_compulsory remember=1\
       && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 4 ;;
      esac
     done
     eval "cc_var_val=\"\$$cc_var_name\""
     case "$cc_is_for" in
     src) tlist_s_cc="" ; tlist_b_cc="" ;;
     utils) tlist_s_cc_native="" ; tlist_b_cc_native="" ;;
     esac
     if [ ! -z "$cc_var_val" ] ; then
      # Test compiler
      pretty_print 1 1 "Testing C compiler:"
      # Check compiler exists
      clear_tags
      CC="$cc_var_val"
      ENVIRONMENT_COMMAND=$user_ENVIRONMENT_COMMAND
      var_CC="$ENVIRONMENT_COMMAND ; $CC"
      set -- $(check_compilers CC)
      (($1==0)) && pretty_print 1 1 "Could not find compiler, try again."\
       && echo && continue
      case "$cc_is_for" in
      src)
       tlist_s_cc="${@:3:$2}" ; tlist_b_cc="${@:$(($2+3))}" ;;
      utils)
       tlist_s_cc_native="${@:3:$2}"
       tlist_b_cc_native="${@:$(($2+3))}" ;;
      esac
      pretty_print 1 3 "- C compiler exists"
      # Test compilation
      run_compiler_test compile cc="$cc_var_val" env="$env_cmd" output=1\
       mpirun="$user_RUN_PARALLEL" use_mpirun="$bin_needs_mpirun"
      t1=$?
      ((t1>1)) && ((ignore_bintest==1)) && t1=0
      if ((t1==1)) ; then
       if [ -z "$env_cmd" ] ; then
        pretty_print 1 1 "This does not appear to be a working C compiler,\
         try again."
       else
        pretty_print 1 1 "This does not appear to be a working C compiler,\
         try again.  Bear in mind that the ENVIRONMENT_COMMAND you have written\
         above might be interfering with the compiler!"
       fi
       echo ; continue
      fi
      pretty_print 1 3 "- C compiler seems to work."
      # Try to see how one should interface C and Fortran
      for c2fflags_var_val in\
       "" "-DF90_CAPITALS"\
       "-DF90_NO_UNDERSCORE" "-DF90_NO_UNDERSCORE -DF90_CAPITALS"\
       "-DF90_DOUBLE_UNDERSCORE" "-DF90_DOUBLE_UNDERSCORE -DF90_CAPITALS" ; do
       eval "$c2fflags_var_name=\"\$c2fflags_var_val\""
       run_compiler_test C_interface f90="$f90_var_val" cc="$cc_var_val"\
        cflags="$c2fflags_var_val" env="$env_cmd" output=1\
        mpirun="$user_RUN_PARALLEL" use_mpirun="$bin_needs_mpirun"
       t1=$?
       ((t1>1)) && ((ignore_bintest==1)) && t1=0
       ((t1==0)) && break
      done
      if ((t1!=0)) ; then
       pretty_print 1 1 "The Fortran and C compilers you have specified do not\
        interface properly. Try again."
       echo
       continue
      fi
      if [ ! -z "$c2fflags_var_val" ] ; then
       pretty_print 1 3 "- C compiler interfaces correctly with Fortran\
        compiler using '$c2fflags_var_val'"
      else
       pretty_print 1 3 "- C compiler interfaces correctly with Fortran\
        compiler without flags"
      fi
      # Try to compile ETIME now, testing for required CFLAGS
      if [ "$cc_is_for" = src ] && [ "$user_NEED_ETIME" = yes ] ; then
       for user_CFLAGS_ETIME in "" "-DCASINO_T3E_MODE" "-DCASINO_NO_ETIME" ; do
        run_compiler_test external_etime f90="$f90_var_val" cc="$cc_var_val"\
         cflags="$c2fflags_var_val $user_CFLAGS_ETIME" env="$env_cmd" output=1\
         mpirun="$user_RUN_PARALLEL" use_mpirun="$bin_needs_mpirun"
        t1=$?
        ((t1>1)) && ((ignore_bintest==1)) && t1=0
        ((t1==0)) && break
       done
       if ((t1!=0)) ; then
        pretty_print 1 1 "The Fortran and C compilers you have specified do not\
         seem to understand each other in producing an ETIME replacement. Try\
         again."
        echo
        continue
       fi
       if [[ "$user_CFLAGS_ETIME" == *NO_ETIME* ]] ; then
        pretty_print 1 3 "- C compiler can only provide a non-functioning ETIME"
       else
        pretty_print 1 3 "- C compiler correctly compiles the ETIME replacement"
       fi
      fi
      # Test SHM support
      if [ "$user_TYPE" != single ] && [ "$cc_is_for" = src ] ; then
       for user_CFLAGS_SHM in "-DSHM_SYSV" "-DSHM_POSIX" ; do
        run_compiler_test shm f90="$f90_var_val" cc="$cc_var_val"\
         fflags="$user_INCLUDE_DIR" lflags="$user_LIB_PATH $user_LDLIBS_all"\
         cflags="$c2fflags_var_val $user_CFLAGS_SHM $user_INCLUDE_DIR"\
         env="$env_cmd" output=1 mpirun="$user_RUN_PARALLEL"
         use_mpirun="$bin_needs_mpirun"
        t1=$?
        ((t1>1)) && ((ignore_bintest==1)) && t1=0
        ((t1==0)) && break
       done
       if ((t1==0)) ; then
        user_SUPPORT_SHM=yes
        pretty_print 1 3 "- C compiler builds the SHM facility with\
         '$user_CFLAGS_SHM'."
       else
        user_SUPPORT_SHM=no
        user_CFLAGS_SHM=""
        pretty_print 1 3 "- C compiler unable to build the SHM facility."
       fi
      fi
      echo
     fi
     break
    done
   done
   # C++ compiler
   step_heading "C++ compiler name"
   while : ; do
    pretty_print 1 1 "A C++ compiler is optional.  Specifying a C++ compiler\
     allows you to build:"
    pretty_print 1 3 "- the optional \"ppconvert\" CASINO-to-ABINIT\
     pseudopotential converter"
    echo
    pretty_print 1 1 "Enter the name of the C++ compiler binary (e.g. 'g++'),\
     or empty to omit:"
    while : ; do
     prompt user_CXX prompt=CXX remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 3 ;;
     esac
    done
    tlist_s_cxx="" ; tlist_b_cxx=""
    if [ ! -z "$user_CXX" ] ; then
     # Test compiler
     pretty_print 1 1 "Testing C++ compiler:"
     if [ "$user_UTILS_MODE" != native ] ; then
      ignore_bintest=$ignore_bintest_requested_utils
     else
      ignore_bintest=$ignore_bintest_requested_src
     fi
     # Check compiler exists
     clear_tags
     CXX=$user_CXX
     ENVIRONMENT_COMMAND=$user_ENVIRONMENT_COMMAND
     var_CXX="$ENVIRONMENT_COMMAND ; $CXX"
     set -- $(check_compilers CXX)
     (($1==0)) && pretty_print 1 1 "Could not find compiler, try again."\
      && echo && continue
     tlist_s_cxx="${@:3:$2}" ; tlist_b_cxx="${@:$(($2+3))}"
     pretty_print 1 3 "- C++ compiler exists"
     # Test compilation
     run_compiler_test compile cxx=$user_CXX env="$user_ENVIRONMENT_COMMAND"\
      output=1
     t1=$?
     ((ignore_bintest==1)) && ((t1>1)) && t1=0
     if ((t1!=0)) ; then
      if [ -z "$user_ENVIRONMENT_COMMAND" ] ; then
       pretty_print 1 1 "This does not appear to be a working C++ compiler,\
        try again."
      else
       pretty_print 1 1 "This does not appear to be a working C++ compiler,\
        try again.  Bear in mind that the ENVIRONMENT_COMMAND you have written\
        above might be interfering with the compiler!"
       echo ; continue
      fi
     fi
     pretty_print 1 3 "- C++ compiler seems to work."
     echo
    fi
    break
   done
   # Autodetect compiler
   step_heading "Pre-set compiler flags and features"
#   tlist_s="$tlist_s_f90 $tlist_s_f90_native $tlist_s_cc $tlist_s_cc_native\
#    $tlist_s_cxx"
#   tlist_b="$tlist_b_f90 $tlist_b_f90_native $tlist_b_cc $tlist_b_cc_native\
#    $tlist_b_cxx"
   pretty_print 1 1 "Auto-detecting compilers:"
   for c in F90 F90_NATIVE CC CC_NATIVE CXX ; do
    [ "$user_UTILS_MODE" != native ] && [[ "$c" == *_NATIVE ]] && continue
    if [ "$user_UTILS_MODE" != native ] ; then
     case "$c" in
     F90)        chuman="Fortran compiler" ; cc=F90 ;;
     F90_NATIVE) continue ;;
     CC)         chuman="C compiler"       ; cc=CC ;;
     CC_NATIVE)  continue ;;
     CXX)        chuman="C++ compiler"     ; cc=CXX ;;
     esac
    else
     case "$c" in
     F90)        chuman="Fortran cross-compiler"  ; cc=F90 ;;
     F90_NATIVE) chuman="Native Fortran compiler" ; cc=F90 ;;
     CC)         chuman="C cross-compiler"        ; cc=CC ;;
     CC_NATIVE)  chuman="Native C compiler"       ; cc=CC ;;
     CXX)        chuman="C++ compiler"            ; cc=CXX ;;
     esac
    fi
    clow=$(uncap $c)
    cclow=$(uncap $cc)
    chosen_compiler=""
    eval "uc=\"\$user_$c\""
    if [ ! -z "$uc" ] ; then
     for file in $cclow/* ; do
      # Load compiler tags.
      comp=${file%.arch}
      comp=${comp##*/}
      echo -n " Probing: $chuman: $comp$el$cr"
      clear_tags
      eval "$cc=\"\$user_$c\""
      case "$c" in
      *_NATIVE) ENVIRONMENT_COMMAND="$user_ENVIRONMENT_COMMAND_NATIVE" ;;
      *) ENVIRONMENT_COMMAND="$user_ENVIRONMENT_COMMAND" ;;
      esac
      eval "var_$cc=\"\$ENVIRONMENT_COMMAND ; \$$cc\""
      eval "tlist_s=\"\$tlist_s_$clow\""
      eval "tlist_b=\"\$tlist_b_$clow\""
      if [ ! -z "$tlist_s" ] ; then
       quickload_tags_scalar "$tlist_s" "$file"
       quickload_tags_block "$tlist_b" "$file"
       set -- $(check_compilers_deep $cc)
       (($1==0)) && continue
       in_line ${cc}_CHECK "$@" || continue # avoid rubbish results
      fi
      chosen_compiler="$comp"
      quickload_tags_scalar "DESCRIPTION" "$file"
      comp_desc="$DESCRIPTION"
      break
     done
     echo -n "$el$cr"
     if [ -z "$chosen_compiler" ] ; then
      pretty_print 1 3 "- $chuman: unknown."
     else
      pretty_print 1 3 "- $chuman: $comp_desc."
     fi
    else
     pretty_print 1 3 "- $cword: not using one."
    fi
    eval "chosen_$clow=\"\$chosen_compiler\""
   done
   echo
   # We force entering compiler options only if we have not chosen an F90
   clear_tags
   if [ ! -z "$chosen_f90" ] ; then
    pretty_print 1 1 "Choose an option:"
    pretty_print 1 5 "[1] Browse the pre-set compiler options, and possibly\
     change them."
    pretty_print 1 5 "[2] Use the pre-set compiler options without reviewing."
    while : ; do
     prompt choice allow="1 2" && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    if ((choice==1)) ; then
     if [ ! -z "$chosen_f90" ] ; then
      quickload_tags_make "F90 FFLAGS_opt FFLAGS_debug FFLAGS0_libs\
       SUPPORT_OPENMP FFLAGS_OPENMP_yes LIB_PATH INCLUDE_DIR LDLIBS_all\
       NEED_ETIME" "f90/$chosen_f90.arch"
      quickload_tags_scalar "F90_VERSION" "f90/$chosen_f90.arch"
      quickload_tags_block "COMMAND_CHECK_F90 COMMAND_CHECK_F90_VERSION"\
       "f90/$chosen_f90.arch"
     fi
     if [ ! -z "$chosen_cc" ] ; then
      quickload_tags_make "CC SUPPORT_SHM CFLAGS_SHM CFLAGS_OPENMP_yes"\
       "cc/$chosen_cc.arch"
      quickload_tags_scalar "CC_VERSION" "cc/$chosen_cc.arch"
      quickload_tags_block "COMMAND_CHECK_CC COMMAND_CHECK_CC_VERSION"\
       "cc/$chosen_cc.arch"
     fi
     if [ ! -z "$chosen_cxx" ] ; then
      quickload_tags_make "CXX" "cxx/$chosen_cxx.arch"
      quickload_tags_scalar "CXX_VERSION" "cxx/$chosen_cxx.arch"
      quickload_tags_block "COMMAND_CHECK_CXX COMMAND_CHECK_CXX_VERSION"\
       "cxx/$chosen_cxx.arch"
     fi
     review_compiler=1
    else
     review_compiler=0
    fi
   else
    review_compiler=1
   fi
   if ((review_compiler>0)) ; then
    # Basic compiler flags
    step_heading "Basic compiler flags"
    if [ "$user_UTILS_MODE" = native ] ; then
     pretty_print 1 1 "NOTE: all options here refer to the cross compilers"
    fi
    pretty_print 1 1 "Enter the F90 compiler flags for full optimization:"
    while : ; do
     prompt user_FFLAGS_opt prompt="FFLAGS(OPT)" default="$FFLAGS_opt"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    pretty_print 1 1 "Enter the F90 compiler flags for debugging:"
    while : ; do
     prompt user_FFLAGS_debug prompt="FFLAGS(DEBUG)" default="$FFLAGS_debug"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    pretty_print 1 1 "Enter the F90 compiler flags for no optimization:"
    [ -z "$FFLAGS0_libs" ] && FFLAGS0_libs="-O0"
    while : ; do
     prompt user_FFLAGS0_libs prompt="FFLAGS(NOOPT)" default="$FFLAGS0_libs"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    # OpenMP
    step_heading "OpenMP configuration"
    suggest_choice=""
    [ "$SUPPORT_OPENMP" = yes ] && suggest_choice=1
    pretty_print 1 1 "The following options are for OpenMP support.  If you\
     don't know what this is, do not know what compiler option is required, or\
     do not wish to use it (the benefits are limited at present), feel free to\
     select option 2."
    echo
    pretty_print 1 1 "Choose an option:"
    pretty_print 1 5 "[1] This F90 compiler supports OpenMP"
    pretty_print 1 5 "[2] This F90 compiler does not support OpenMP"
    while : ; do
     prompt choice default="$suggest_choice" allow="1 2" && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    case "$choice" in
    1)
     user_SUPPORT_OPENMP=yes
     pretty_print 1 1 "Enter the Fortran compiler flags required to compile\
      OpenMP code (e.g., '-mp', '-openmp', '-fopenmp -lsci_quadcore_mp', ...):"
     while : ; do
      prompt user_FFLAGS_OPENMP_yes prompt=FFLAGS_OPENMP\
       default="$FFLAGS_OPENMP_yes" remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
     if [ "$user_SUPPORT_SHM" = yes ] ; then
      pretty_print 1 1 "Enter the C compiler flags required to compile\
       OpenMP code (e.g., '-mp', '-openmp', '-fopenmp -lsci_quadcore_mp', ...):"
      while : ; do
       prompt user_CFLAGS_OPENMP_yes prompt=CFLAGS_OPENMP\
        default="$CFLAGS_OPENMP_yes" remember=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 2 ;;
       esac
      done
     fi ;;
    2)
     user_SUPPORT_OPENMP=no
     user_FFLAGS_OPENMP_yes="" ; user_CFLAGS_OPENMP_yes="" ;;
    esac
   fi

   case "$user_TYPE" in
   single) next_stage=ARCH_DETAILS ;;
   *)      next_stage=TOPOLOGY ;;
   esac ;;

  TOPOLOGY)
   # Gather TYPE-specific information
   case "$user_TYPE" in
   single) # should not be run
    pretty_print 1 1 "No configuration required for single-core machine." ;;
   parallel)
    step_heading "Running parallel jobs"
    # CORES_PER_NODE
    suggest_CPN=0
    if [ ! -z "$chosen_os" ] ; then
     clear_tags
     load_tags "os/$chosen_os.arch"
     eval_tag CORES_PER_NODE
     check_number_N "$CORES_PER_NODE" && suggest_CPN=$CORES_PER_NODE
    fi
    if ((suggest_CPN>0)) ; then
     pretty_print 1 1 "Pre-set OS definitions contain auto-detection code that\
      can determine the number of cores ($suggest_CPN), no set-up required."
     echo
    else
     pretty_print 1 1 "No OS autodetection.  How many cores does this machine\
      have?"
     while : ; do
      prompt user_CORES_PER_NODE prompt=CORES_PER_NODE force_integer=1\
       remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
    fi
    # RUN_PARALLEL
    pretty_print 1 1 "Enter the command required to run an MPI calculation.\
     Here you must use the following variables:"
    pretty_print 1 3 "- &NPROC& expands to the number of MPI processes"
    pretty_print 1 3 "- &BINARY& expands to the name of the binary"
    while : ; do
     while : ; do
      prompt user_RUN_PARALLEL prompt=RUN_PARALLEL\
       default="mpirun -np &NPROC& &BINARY&" non_empty=1 remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 3 ;;
      esac
     done
     ( [[ "$user_RUN_PARALLEL" != *"&NPROC&"* ]]\
      || [[ "$user_RUN_PARALLEL" != *"&BINARY&"* ]] )\
      && pretty_print 1 1 "The command must contain the &NPROC& and &BINARY&\
      variables, try again." && echo && continue
     break
    done ;;
   cluster)
    step_heading "Machine topology"
    # CORES_PER_NODE
    suggest_CPN=0
    if [ ! -z "$chosen_os" ] ; then
     clear_tags
     load_tags "os/$chosen_os.arch"
     eval_tag CORES_PER_NODE
     check_number_N "$CORES_PER_NODE" && suggest_CPN=$CORES_PER_NODE
    fi
    if ((suggest_CPN>0)) ; then
     pretty_print 1 1 "Pre-set OS definitions contain auto-detection code that\
      can determine the number of cores in log-in nodes ($suggest_CPN), no\
      set-up required."
     echo
    else
     pretty_print 1 1 "No OS autodetection.  How many cores are there in the\
      log-in node of this machine?"
     while : ; do
      prompt user_CORES_PER_NODE prompt=CORES_PER_NODE force_integer=1\
       remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
    fi
    pretty_print 1 1 "Computational nodes may differ from log-in nodes in some\
     machines, although they are often identical.  How many cores are there in\
     each computational node of this machine?"
    if ((suggest_CPN>0)) ; then
     while : ; do
      prompt user_CORES_PER_NODE_CLUSTER prompt=CORES_PER_NODE_CLUSTER\
       force_integer=1 default="$suggest_CPN" remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
    else
     while : ; do
      prompt user_CORES_PER_NODE_CLUSTER prompt=CORES_PER_NODE_CLUSTER\
       force_integer=1 default="$CORES_PER_NODE" remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
    fi
    step_heading "Running parallel jobs on the log-in node"
    pretty_print 1 1 "Enter the command required to run an MPI calculation on\
     the log-in node. Here you must use the following variables:"
    pretty_print 1 3 "- &NPROC& expands to the number of MPI processes"
    pretty_print 1 3 "- &BINARY& expands to the name of the binary"
    while : ; do
     while : ; do
      prompt user_RUN_PARALLEL prompt=RUN_PARALLEL\
       default="mpirun -np &NPROC& &BINARY&" non_empty=1 remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 3 ;;
      esac
     done
     ( [[ "$user_RUN_PARALLEL" != *"&NPROC&"* ]]\
      || [[ "$user_RUN_PARALLEL" != *"&BINARY&"* ]] )\
      && pretty_print 1 1 "The command must contain the &NPROC& and &BINARY&\
      variables, try again." && echo && continue
     break
    done
    step_heading "Pre-set queueing system configuration"
    pretty_print 1 1 "Choose an option."
    allow=""
    i=0 ; for file in queue/* ; do i=$((i+1))
     clear_tags
     quickload_tags_scalar "DESCRIPTION QUEUEING_SYSTEM" "$file"
     pretty_print 1 5 "[$i] The queueing system on this machine is\
      $DESCRIPTION"
     file="${file#*/}" ; file="${file%.arch}"
     choice_queue[$i]="$file"
     choice_QUEUING_SYSTEM[$i]="$QUEUEING_SYSTEM"
     allow="$allow $i"
    done
    i=$((i+1)) ; allow="$allow $i"
    pretty_print 1 5 "[$i] The queueing system on this machine is none of the\
     above"
    while : ; do
     prompt choice allow="$allow" && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
    if ((choice<i)) ; then
     chosen_queue="${choice_queue[$choice]}"
     user_QUEUEING_SYSTEM="${choice_QUEUEING_SYSTEM[$choice]}"
    elif ((choice==i-1)) ; then
     chosen_queue=""
    fi
    if [ ! -z "$chosen_queue" ] ; then
     pretty_print 1 1 "Choose an option:"
     pretty_print 1 5 "[1] Browse the pre-set queueing-system options, and\
      possibly change them."
     pretty_print 1 5 "[2] Use the pre-set queueing-system options without\
      reviewing."
     while : ; do
      prompt choice allow="1 2" && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
     case "$choice" in
     1) review_queue=1 ;;
     2) review_queue=0 ;;
     esac
    else
     step_heading "New queueing system configuration"
     pretty_print 1 1 "Enter the name of this queueing system at the prompt\
      below."
     while : ; do
      prompt user_QUEUEING_SYSTEM prompt=QUEUEING_SYSTEM non_empty=1\
       remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
     review_queue=1
    fi
    if ((review_queue>0)) ; then
     step_heading "Choose an editor"
     # Prompt for editor
     pretty_print 1 1 "Clusters usually require a batch-submission script in\
      order to run jobs. Now you will have to write the template batch\
      submission script that will be used for all CASINO jobs on this machine.\
      This will be done in an external editor."
     echo
     pretty_print 1 1 "Enter the executable name of your preferred editor at\
      the prompt below:"
     while : ; do
      while : ; do
       prompt pref_editor non_empty=1 default="$EDITOR" remember=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      ! type -P "$pref_editor" >& /dev/null && pretty_print 1 1 "$pref_editor\
       not found, try again." && echo && continue
      break
     done
     cut_here_string="--- WRITE BELOW THIS LINE ---"
     # Prompt for head of submission script
     step_heading "Header of submission script"
     tmp_dir="$(cd "$($mktemp_d_command)" >& /dev/null ; pwd)"
     bscript="$tmp_dir/BATCH_SCRIPT_HEAD"
     touch "$bscript"
     cat >> "$bscript" <<___EOF
 Write the HEADER of the submission script below, NOT including the command
 required to run the calculation (typically mpirun or aprun); we will do that
 a bit later.  This typically consists of queueing-system directives, e.g.,
 a PBS directive looks like this:
  #PBS -l walltime=&WALLTIME&:nodes=&NNODE&:ppn=&PPN&
 while an SGE directive looks like this:
  #\$-N &SCRIPT&

 The best thing to do here is to copy the header of a sample submission script
 from the cluster's documentation (there might be a 'Running Jobs' section in
 the cluster's website, for example) and paste it below, and then replace the
 number of processes, walltime, etc, with the variables listed below.  Once
 you are done, save the file and exit the editor, and the configuration script
 will resume.

 You have the following variables available; use them as required:
 - &OUT&      : output file; redirect standard output and standard error here
 - &NPROC&    : number of MPI processes to run
 - &NNODE&    : number of (multi-core) nodes to use
 - &NCORE&    : number of cores to use
 - &PPN&      : MPI processes per node
 - &TPP&      : OpenMP threads per MPI process
 - &TPN&      : OpenMP threads per node
 - &WALLTIME& : requested walltime, in any format (this will be defined later)
 - &BINARY&   : binary to run
 - &SCRIPT&   : name of the script - use this to name the job

 Other variables:
 - &ENV.VARIABLE&
     Use this if you need to access the value of an environment variable in the
     context of the log-in node.  E.g., &ENV.PATH& would expand to the \$PATH
     environment variable in the context of the log-in node.
 - &USER.VARIABLE&
     Use this to define custom variables which will become available for users
     to set on the command line of the runscript.  Common examples are
     &USER.QUEUE&, &USER.ACCOUNT& and &USER.MEM&.  NB, &USER.QUEUE& is
     treated specially by this script.

___EOF
     if [ ! -z "$user_SCRIPT_HEAD" ] ; then
      cat >> "$bscript" <<___EOF
 The previously-entered submission script header has been pre-loaded below.

___EOF
     elif [ ! -z "$chosen_queue" ] ; then
      cat >> "$bscript" <<___EOF
 The default submission script header for the queueing system you have chosen
 has been pre-loaded below.

___EOF
     fi
     echo "$cut_here_string" >> "$bscript"
     clear_tags
     if [ ! -z "$user_SCRIPT_HEAD" ] ; then
      print_block_tag user_SCRIPT_HEAD >> "$bscript"
     elif [ ! -z "$chosen_queue" ] ; then
      quickload_tags_block SCRIPT_HEAD "queue/$chosen_queue.arch"
      [ -z "$SCRIPT_HEAD" ] || print_block_tag SCRIPT_HEAD >> "$bscript"
     fi
     echo -n " Opening editor to edit batch script header"
     echo -n . ; sleep 1 ; echo -n . ; sleep 1 ; echo -n . ; sleep 1 ; echo
     $pref_editor "$bscript"
     user_SCRIPT_HEAD[0]=0
     skip=1
     {
      while : ; do
       IFS=""
       read -r line || break
       IFS="$IFS_save"
       if ((skip==1)) ; then
        [ "$line" = "$cut_here_string" ] && skip=0
       else
        user_SCRIPT_HEAD[0]=$((${user_SCRIPT_HEAD[0]}+1))
        user_SCRIPT_HEAD[${user_SCRIPT_HEAD[0]}]="$line"
       fi
      done
     } < "$bscript"
     IFS="$IFS_save"
     rm -rf "$tmp_dir"
     pretty_print 1 1 "Loaded SCRIPT_HEAD." ; echo
     # Prompt for run command of submission script
     step_heading "Run section of submission script"
     tmp_dir="$(cd "$($mktemp_d_command)" >& /dev/null ; pwd)"
     bscript="$tmp_dir/BATCH_SCRIPT_RUN"
     touch "$bscript"
     cat >> "$bscript" <<___EOF
 Write the run command to use in the submission script below, NOT including the
 header which you have already written.  The run command should consist of a
 single line with mpirun, aprun or the machine's equivalent.

 Again, copy the section of a sample submission script from the cluster's
 documentation that executes this command and paste it below, then replace the
 number of processes, name of the binary, etc, with the variables listed below.
 Once you are done, save the file and exit the editor, and the configuration
 script will resume.

 You have the following variables available; use them as required:
 - &OUT&      : output file; redirect standard output and standard error here
 - &NPROC&    : number of MPI processes to run
 - &NNODE&    : number of (multi-core) nodes to use
 - &NCORE&    : number of cores to use
 - &PPN&      : MPI processes per node
 - &TPP&      : OpenMP threads per MPI process
 - &TPN&      : OpenMP threads per node
 - &WALLTIME& : requested walltime, in any format (this will be defined later)
 - &BINARY&   : binary to run
 - &SCRIPT&   : name of the script - use this to name the job

 Other variables:
 - &ENV.VARIABLE&
     Use this if you need to access the value of an environment variable in the
     context of the log-in node.  E.g., &ENV.PATH& would expand to the \$PATH
     environment variable in the context of the log-in node.
 - &USER.VARIABLE&
     Use this to define custom variables which will become available for users
     to set on the command line of the runscript.  Common examples are
     &USER.QUEUE&, &USER.ACCOUNT& and &USER.MEM&.  NB, &USER.QUEUE& is
     treated specially by this script.

___EOF
     if [ ! -z "$user_SCRIPT_RUN" ] ; then
      cat >> "$bscript" <<___EOF
 The previously-entered submission script header has been pre-loaded below.

___EOF
     elif [ ! -z "$chosen_queue" ] ; then
      cat >> "$bscript" <<___EOF
 The default submission run command for the queueing system you have chosen
 has been pre-loaded below.

___EOF
     fi
     echo "$cut_here_string" >> "$bscript"
     clear_tags
     if [ ! -z "$user_SCRIPT_RUN" ] ; then
      print_block_tag user_SCRIPT_RUN >> "$bscript"
     elif [ ! -z "$chosen_queue" ] ; then
      quickload_tags_block SCRIPT_RUN "queue/$chosen_queue.arch"
      [ -z "$SCRIPT_RUN" ] || print_block_tag SCRIPT_RUN >> "$bscript"
     fi
     echo -n " Opening editor to edit batch script run command"
     echo -n . ; sleep 1 ; echo -n . ; sleep 1 ; echo -n . ; sleep 1 ; echo
     $pref_editor "$bscript"
     user_SCRIPT_RUN[0]=0
     skip=1
     {
      while : ; do
       IFS=""
       read -r line || break
       IFS="$IFS_save"
       if ((skip==1)) ; then
        [ "$line" = "$cut_here_string" ] && skip=0
       else
        user_SCRIPT_RUN[0]=$((${user_SCRIPT_RUN[0]}+1))
        user_SCRIPT_RUN[${user_SCRIPT_RUN[0]}]="$line"
       fi
      done
     } < "$bscript"
     IFS="$IFS_save"
     rm -rf "$tmp_dir"
     pretty_print 1 1 "Loaded SCRIPT_RUN." ; echo
     # Prompt for submission command.
     step_heading "Submission command"
     pretty_print 1 1 "Batch scripts are usually submitted to the queue using\
      a command named qsub, bsub, llsubmit, or similar.  Typical submission\
      commands are:"
     pretty_print 1 3 "- qsub &SCRIPT&"
     pretty_print 1 3 "- bsub < &SCRIPT&"
     pretty_print 1 3 "- llsubmit &SCRIPT&"
     pretty_print 1 1 "You will need to provide the appropriate command below."
     echo
     pretty_print 1 1 "You have the following variables available; use them as\
      required (you should only need &SCRIPT&, though):"
     pretty_print 1 3 "- &OUT&      : output file; redirect standard output\
      and standard error here"
     pretty_print 1 3 "- &NPROC&    : number of MPI processes to run"
     pretty_print 1 3 "- &NNODE&    : number of (multi-core) nodes to use"
     pretty_print 1 3 "- &NCORE&    : number of cores to use"
     pretty_print 1 3 "- &PPN&      : MPI processes per node"
     pretty_print 1 3 "- &TPP&      : OpenMP threads per MPI process"
     pretty_print 1 3 "- &TPN&      : OpenMP threads per node"
     pretty_print 1 3 "- &WALLTIME& : requested walltime, in any format (this\
      will be defined later)"
     pretty_print 1 3 "- &BINARY&   : binary to run"
     pretty_print 1 3 "- &SCRIPT&   : name of the script"
     echo
     pretty_print 1 1 "Other variables:"
     pretty_print 1 3 "- &ENV.VARIABLE&"
     pretty_print 5 5 "Use this if you need to access the value of an\
      environment variable in the context of the log-in node.  E.g.,\
      &ENV.PATH& would expand to the \$PATH environment variable in the\
      context of the log-in node."
     pretty_print 1 3 "- &USER.VARIABLE&"
     pretty_print 5 5 "Use this to define custom variables which will become\
      available for users to set on the command line of the runscript.  Common\
      examples are &USER.QUEUE&, &USER.ACCOUNT& and &USER.MEM&.  NB,\
      &USER.QUEUE& is treated specially by this script."
     echo
     pretty_print 1 1 "Enter the command required to submit a batch script\
      named &SCRIPT&:"
     clear_tags
     [ -z "$chosen_queue" ] || quickload_tags_scalar SUBMIT_SCRIPT\
      "queue/$chosen_queue.arch"
     while : ; do
      prompt user_SUBMIT_SCRIPT prompt=SUBMIT_SCRIPT default="$SUBMIT_SCRIPT"\
       remember=1 non_empty=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
    fi
    step_heading "Machine policies"
    if ((review_queue==0)) ; then # implies [ ! -z "$chosen_queue" ]
     # Find which variables have been used
     clear_tags
     quickload_tags_block "SCRIPT_HEAD SCRIPT_RUN" "queue/$chosen_queue.arch"
     quickload_tags_scalar SUBMIT_SCRIPT "queue/$chosen_queue.arch"
     cluster_variables="$SCRIPT_HEAD_deps $SCRIPT_RUN_deps $SUBMIT_SCRIPT_deps"
    else
     # Find which variables have been used
     cluster_variables=""
     if [ ! -z "$user_SCRIPT_HEAD" ] ; then
      i=0 ; while ((i<${user_SCRIPT_HEAD[0]})) ; do i=$((i+1))
       cluster_variables="$cluster_variables $(variables_in_line\
        "${user_SCRIPT_HEAD[$i]}")"
      done
     fi
     if [ ! -z "$user_SCRIPT_RUN" ] ; then
      i=0 ; while ((i<${user_SCRIPT_RUN[0]})) ; do i=$((i+1))
       cluster_variables="$cluster_variables $(variables_in_line\
        "${user_SCRIPT_RUN[$i]}")"
      done
     fi
     cluster_variables="$cluster_variables $(variables_in_line\
      "$user_SUBMIT_SCRIPT")"
    fi
    cluster_variables="$(uniq_list $cluster_variables)"
    # User variables
    cluster_variables_user=""
    cluster_variables_internal=""
    for var in $cluster_variables ; do
     case "$var" in
     USER.*)
      var="${var#*.}"
      cluster_variables_user="$cluster_variables_user $var"
      pretty_print 1 1 "You have specified the &USER.$var& variable.  This\
       variable will be available on the command line when calculations are\
       run.  Choose an option:"
      pretty_print 1 5 "[1] Specify a strict list of possible values for\
       &USER.$var&, the first of which will be the default"
      pretty_print 1 5 "[2] Specify just a default value &USER.$var&, allowing\
       any other values to be used"
      pretty_print 1 5 "[3] Specify neither, requiring the value to be\
       specified every time a calculation is run, and allowing any value to be\
       used"
      while : ; do
       prompt choice allow="1 2 3" && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      case "$choice" in
      1)
       local user_USER_ALLOWED_$var
       pretty_print 1 1 "Enter the possible values of &USER.$var& as a\
        space-separated list, the first of which will be the default"
       while : ; do
        prompt user_USER_ALLOWED_$var prompt="USER.ALLOWED.$var" non_empty=1\
         remember=1 && break
        meta_menu
        case $? in
        0) : ;;
        1) return 1 ;;
        2) continue 3 ;;
        esac
       done ;;
      2)
       local user_USER_DEFAULT_$var
       pretty_print 1 1 "Enter the default values of &USER.$var&"
       while : ; do
        prompt user_USER_DEFAULT_$var prompt="USER.DEFAULT.$var" non_empty=1\
         remember=1 && break
        meta_menu
        case $? in
        0) : ;;
        1) return 1 ;;
        2) continue 3 ;;
        esac
       done ;;
      esac
      suggest_desc=""
      case "$var" in
      QUEUE) suggest_desc="Name of queue under which to run the job" ;;
      ACCOUNT) suggest_desc="Name of account under which to run the job" ;;
      MEM) suggest_desc="Requested memory per node" ;;
      esac
      local user_USER_DESCRIPTION_$var
      pretty_print 1 1 "Enter a description for &USER.$var&, which will be\
       displayed when the --help option is passed to the runscript."
      while : ; do
       prompt user_USER_DESCRIPTION_$var prompt="USER.DESCRIPTION.$var"\
        default="$suggest_desc" non_empty=1 remember=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done ;;
     INTERNAL.*)
      var="${var#*.}"
      cluster_variables_internal="$cluster_variables_internal $var" ;;
     ENV.*)
      ;;
     *.*) # this is an error - should do something
      ;;
     *) # should check if this is not a defined variable and do something
      ;;
     esac
    done
    # Size and time limits
    if in_line WALLTIME $cluster_variables ; then
     pretty_print 1 1 "What is the format in which walltimes should be\
      specified?  Use D, H, M, S to denote days, hours, minutes and seconds,\
      and DD, DDD etc if you need zero-padding to 2-, 3- etc digits."
     echo
     pretty_print 1 1 "Example: a walltime of 5h43m31s would see the following\
      expansions for different formats:"
     pretty_print 1 3 "- H:MM:SS -> 5:43:31"
     pretty_print 1 3 "- DD:HH:M:SSSS -> 00:05:43:0031"
     pretty_print 1 3 "- MMMM minutes SS seconds -> 0343 minutes 31 seconds"
     echo
     pretty_print 1 1 "Enter the time format string:"
     while : ; do
      prompt user_TIME_FORMAT prompt=TIME_FORMAT default="H:MM:SS"\
       non_empty=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
     pretty_print 1 1 "The CASINO runscript needs to know the limits on\
      walltime and number of cores for jobs run on this machine.  There are\
      different types of job-size policies on different machines; choose the\
      option below that better describes this machine's policy."
     echo
     if [ ! -z "$user_USER_ALLOWED_QUEUE" ] ; then
      pretty_print 1 1 "Choose an option:"
      pretty_print 1 5 "[1] This machine has a number-of-cores limit only, and\
       no walltime limit"
      pretty_print 1 5 "[2] This machine has one walltime limit and one\
       number-of-cores limit"
      pretty_print 1 5 "[3] This machine has one number-of-cores limit for\
       each queue, and no walltime limit"
      pretty_print 1 5 "[4] This machine has one walltime limit and one\
       number-of-cores limit for each queue"
      pretty_print 1 5 "[5] This machine has walltime limits depending on the\
       number of cores used"
      while : ; do
       prompt choice allow="1 2 3 4 5" && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 2 ;;
       esac
      done
      case "$choice" in
      1) limits_type=ncore ;;
      2) limits_type=indep ;;
      3) limits_type=queue_ncore ;;
      4) limits_type=queue ;;
      5) limits_type=interdep ;;
      esac
     else
      pretty_print 1 1 "Choose an option:"
      pretty_print 1 5 "[1] This machine has a number-of-cores limit only, and\
       no walltime limit"
      pretty_print 1 5 "[2] This machine has one walltime limit and one\
       number-of-cores limit"
      pretty_print 1 5 "[3] This machine has walltime limits depending on the\
       number of cores used"
      while : ; do
       prompt choice allow="1 2 3" && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 2 ;;
       esac
      done
      case "$choice" in
      1) limits_type=ncore ;;
      2) limits_type=indep ;;
      3) limits_type=interdep ;;
      esac
     fi
    else
     pretty_print 1 1 "You have not used variable &WALLTIME&, so CASINO does\
      not need to know whether there are any walltime limits in place on this\
      machine."
     echo
     if [ ! -z "$user_USER_ALLOWED_QUEUE" ] ; then
      pretty_print 1 1 "Choose an option:"
      pretty_print 1 5 "[1] This machine has a global number-of-cores limit"
      pretty_print 1 5 "[2] This machine has one number-of-cores limit for\
       each queue"
      while : ; do
       prompt choice allow="1 2" && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 2 ;;
       esac
      done
      case "$choice" in
      1) limits_type=ncore ;;
      2) limits_type=queue_ncore ;;
      esac
     else
      limits_type=ncore
     fi
    fi
    step_heading "Machine job limits"
    case "$limits_type" in
    ncore)
     pretty_print 1 1 "Enter the maximum number of cores a single job can use\
      on this machine:"
     while : ; do
      prompt user_MAX_NCORE prompt=MAX_NCORE force_integer=1 remember=1\
       && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done ;;
    indep)
     pretty_print 1 1 "Enter the maximum number of cores a single job can use\
      on this machine:"
     while : ; do
      prompt user_MAX_NCORE prompt=MAX_NCORE force_integer=1\
       remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done
     pretty_print 1 1 "Enter the maximum walltime a single job can use on this\
      machine (use the format <days>d<hours>h<minutes>m<seconds>s, e.g.\
      1d7h5s, 36h1m):"
     while : ; do
      prompt user_MAX_WALLTIME prompt=MAX_WALLTIME non_empty=1\
       remember=1 && break
      meta_menu
      case $? in
      0) : ;;
      1) return 1 ;;
      2) continue 2 ;;
      esac
     done ;;
    queue_ncore)
     user_MAX_NCORE_per_queue=""
     for queue in $user_USER_ALLOWED_QUEUE ; do
      pretty_print 1 1 "Enter the number-of-cores limit for queue '$queue':"
      while : ; do
       prompt choice prompt="MAX_NCORE($queue)" force_integer=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      user_MAX_NCORE_per_queue="$user_MAX_NCORE_per_queue $choice"
     done ;;
    queue)
     user_MAX_NCORE_per_queue=""
     user_MAX_WALLTIME_per_queue=""
     for queue in $user_USER_ALLOWED_QUEUE ; do
      pretty_print 1 1 "Enter the number-of-cores limit for queue '$queue':"
      while : ; do
       prompt choice prompt="MAX_NCORE($queue)" force_integer=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      user_MAX_NCORE_per_queue="$user_MAX_NCORE_per_queue $choice"
      pretty_print 1 1 "Enter the walltime limit for queue '$queue':"
      while : ; do
       prompt choice prompt="MAX_WALLTIME($queue)" non_empty=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      user_MAX_WALLTIME_per_queue="$user_MAX_WALLTIME_per_queue $choice"
     done ;;
    interdep)
     pretty_print 1 1 "In this set-up there are 'bands' of number of cores,\
      each of which has a different walltime limit.  You will need to give the\
      maximum number of cores in each band (in increasing order), followed by\
      the walltime limit of that band (in the format\
      <days>d<hours>h<minutes>m<seconds>s, e.g. 1d7h5s, 36h1m).  Enter 0 as\
      the number of cores to finish setting up bands."
     echo
     user_MAX_NCORE_interdep=""
     user_MAX_WALLTIME_interdep=""
     last_ncore=0
     i=0 ; while : ; do i=$((i+1))
      pretty_print 1 1 "Enter the maximum number of cores in band #$i (0 to\
       finish):"
      while : ; do
       prompt choice prompt="MAX_NCORE($i)" force_integer=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      ((choice==0)) && break
      ((last_ncore>choice)) && pretty_print 1 1 "The bands must be given in\
       increasing order, try again." && echo && i=$((i-1)) && continue
      user_MAX_NCORE_interdep="$user_MAX_NCORE_interdep $choice"
      pretty_print 1 1 "Enter the walltime limit for band #$i:"
      while : ; do
       prompt choice prompt="MAX_WALLTIME($i)" non_empty=1 && break
       meta_menu
       case $? in
       0) : ;;
       1) return 1 ;;
       2) continue 3 ;;
       esac
      done
      user_MAX_WALLTIME_interdep="$user_MAX_WALLTIME_interdep $choice"
     done ;;
    esac ;;
   esac

   next_stage=ARCH_DETAILS ;;

  ARCH_DETAILS)
   # Prepare suggested CASINO_ARCH name
   suggest_system="$(uncap $host_KERNEL)"
   suggest_system="${suggest_system//\//_}"
   case "$suggest_system" in
   linux) suggest_system=linuxpc ;;
   darwin) suggest_system=macos ;;
   sunos) suggest_system=sun ;;
   osf1) suggest_system=alpha ;;
   esac
   if [ -z "$chosen_f90" ] ; then
    suggest_compiler="$user_F90"
   else
    suggest_compiler="$chosen_f90"
   fi
   case "$user_TYPE" in
   single)
    arch_form="<system>-<compiler>"
    suggest_arch="$suggest_system-$suggest_compiler" ;;
   parallel)
    arch_form="<system>-<compiler>-parallel"
    suggest_arch="$suggest_system-$suggest_compiler-parallel" ;;
   cluster)
    arch_form="<system>-<compiler>-<queueing-system>-parallel"
    suggest_queue="$chosen_queue"
    [ -z "$suggest_queue" ] && suggest_queue="$(uncap $user_QUEUEING_SYSTEM)"
    [ -z "$suggest_queue" ] && suggest_queue=cluster
    suggest_arch="$suggest_system-$suggest_compiler-$suggest_queue-parallel" ;;
   esac
   # Ask whether this is host-specific, ask for CASINO_ARCH and description
   step_heading "Hostnames"
   pretty_print 1 1 "You can specify a set of machine hostnames if this set-up\
    is specific to a small number of machines.  By doing so you allow other\
    potential users of CASINO on those machines to quickly find this\
    CASINO_ARCH."
   echo
   pretty_print 1 1 "You can use bash-style patterns; the following are\
    particularly useful:"
   pretty_print 1 3 "- * will match zero or more characters"
   pretty_print 1 3 "- ? will match zero or one character"
   echo
   pretty_print 1 1 "Examples of what to write here:"
   pretty_print 1 3 "- pinta??"
   pretty_print 1 3 "- bindloe??, *.hpc.cam.ac.uk"
   echo
   pretty_print 1 1 "Enter the comma-separated list of hostnames or hostname\
    patterns that this CASINO_ARCH is applicable to, or empty to omit:"
   while : ; do
    while : ; do
     prompt user_HOSTNAME prompt=HOSTNAME default="$host_HOSTNAME"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 3 ;;
     esac
    done
    if [ ! -z "$user_HOSTNAME" ] ; then
     clear_tags
     HOSTNAME="$user_HOSTNAME"
     set -- $(check_host)
     if [ "$1" = 0 ] ; then
      pretty_print 1 1 "Current hostname ($host_HOSTNAME) does not match, try\
       again."
      echo ; continue
     fi
    fi
    break
   done
   if [ ! -z "$user_HOSTNAME" ] ; then
    arch_form="$arch_form.<host>"
    suggest_arch="$suggest_arch.${host_HOSTNAME%%.*}"
   fi
   step_heading "CASINO_ARCH name"
   pretty_print 1 1 "Now we need the name for this CASINO_ARCH, which should\
    be roughly of the form $arch_form."
   echo
   pretty_print 1 1 "Enter the name of this CASINO_ARCH:"
   while : ; do
    while : ; do
     prompt var=user_CASINO_ARCH prompt=CASINO_ARCH default="$suggest_arch"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 3 ;;
     esac
    done
    [ -z "$user_CASINO_ARCH" ] && user_CASINO_ARCH="$suggest_arch"
    [ -z "$user_CASINO_ARCH" ] && pretty_print 1 1 "Must provide a name, try\
     again." && echo && continue
    [[ "$user_CASINO_ARCH" != +([[:alnum:]-._]) ]]\
     && pretty_print 1 1 "The name may only contain letters, numbers, '-',\
     '.', or '_', try again." && echo && continue
    [ -e "$user_CASINO_ARCH.arch" ]\
     && pretty_print 1 1 "This name is already in use, try a different one."\
     && echo && continue
    break
   done
   step_heading "CASINO_ARCH description"
   pretty_print 1 1 "Enter a one-line description of this CASINO_ARCH, for\
    documentative purposes:"
   while : ; do
    prompt user_DESCRIPTION prompt=DESCRIPTION non_empty=1 remember=1 && break
    meta_menu
    case $? in
    0) : ;;
    1) return 1 ;;
    2) continue 2 ;;
    esac
   done
   # Ask user to become maintainer
   step_heading "Your details"
   pretty_print 1 1 "You should write your details so that people can contact\
    you regarding this machine and its configuration."
   echo
   pretty_print 1 1 "Enter your name in the following prompt, and your e-mail\
    address in the next:"
   if type -P getent >& /dev/null ; then
    while : ; do
     prompt name prompt=NAME\
      default="$(getent passwd $USER | cut -d ':' -f 5 | cut -d ',' -f 1)"\
      remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
   else
    while : ; do
     prompt name prompt=NAME remember=1 && break
     meta_menu
     case $? in
     0) : ;;
     1) return 1 ;;
     2) continue 2 ;;
     esac
    done
   fi
   while : ; do
    prompt email prompt=EMAIL default="$USER@$host_HOSTNAME" remember=1\
     && break
    meta_menu
    case $? in
    0) : ;;
    1) return 1 ;;
    2) continue 2 ;;
    esac
   done
   if [ ! -z "$name" ] && [ ! -z "$email" ] ; then
    user_MAINTAINER="$name <$email>"
   elif [ ! -z "$name" ] ; then
    user_MAINTAINER="$name"
   elif [ ! -z "$email" ] ; then
    user_MAINTAINER="<$email>"
   else
    user_MAINTAINER="unmaintained"
    pretty_print 1 1 "No details, resulting CASINO_ARCH will be flagged as\
     unmaintained."
    echo
   fi

   next_stage=WRITE ;;

  WRITE)
   break ;;

  *)
   pretty_print 1 1 "Reached unknown stage '$current_stage', aborting."
   return 1 ;;

  esac

  # Keep history of previous stages so we can go back using meta-menu.
  stages_done="$stages_done $current_stage"

 done

 # Write .arch file.
 step_heading "Writing CASINO_ARCH [automatic]"
 echo -n " Writing CASINO/arch/data/$user_CASINO_ARCH.arch...$el$cr"
 arch_file=$user_CASINO_ARCH.arch
 touch "$arch_file"
 cmt1=0
 wrt() {
  ((cmt1==1)) && echo >> "$arch_file" && echo "# $cmt" >> "$arch_file" ; cmt1=0
  echo "$1" >> "$arch_file"
 }
 wrt_host() {
  local varname="$1" host_val include_val
  eval "host_val=\"\$host_$varname\" ; include_val=\"\$$varname\""
  [ "$host_val" = "$include_val" ] && return 0
  wrt "#-! $varname: $host_val"
 }
 wrt_scalar() {
  local varname="$1" user_val
  eval "user_val=\"\$user_$varname\""
  wrt "#-! $varname: $user_val"
 }
 wrt_make() {
  local varname="$1" valdefault="$2" user_val include_val
  eval "user_val=\"\$user_$varname\" ; include_val=\"\$$varname\""
  [ "$user_val" = "$include_val" ] && return 0
  if [ ! -z "$valdefault" ] ; then
   [ "$user_val" = "$valdefault" ] && [ -z "$include_val" ] && return 0
   [ -z "$user_val" ] && [ "$include_val" = "$valdefault" ] && return 0
   [ -z "$user_val" ] && user_val="$valdefault"
  fi
  wrt "$varname = $user_val"
 }
 wrt_make_native() {
  local varname="$1" valdefault="$2" user_val include_val
  eval "user_val=\"\$user_${varname}_NATIVE\" ; include_val=\"\$$varname\""
  [ -z "$user_val" ] && user_val="$include_val"
  [ "$user_val" = "$valdefault" ] && return 0
  [ -z "$user_val" ] && [ ! -z "$valdefault" ] && return 0
  wrt "${varname}_NATIVE = $user_val"
 }
 wrt_scalar DESCRIPTION
 wrt_scalar MAINTAINER
 user_DATE="$(date +%m.%Y)"
 wrt_scalar DATE
 clear_tags
 [ -z "$chosen_os" ] || load_tags "os/$chosen_os.arch"
 wrt_host OS
 wrt_host KERNEL
 wrt_host ARCH
 [ -z "$user_HOSTNAME" ] || wrt_scalar HOSTNAME
 case "$user_TYPE" in
 parallel)
  [ -z "$user_CORES_PER_NODE" ] || wrt_scalar CORES_PER_NODE
  [ ! -z "$user_RUN_PARALLEL" ]\
   && [ "$user_RUN_PARALLEL" != "mpirun -np &NPROC& &BINARY&" ]\
   && wrt_scalar RUN_PARALLEL ;;
 cluster)
  [ -z "$chosen_queue" ] && wrt_scalar QUEUEING_SYSTEM
  [ -z "$user_CORES_PER_NODE" ] || wrt_scalar CORES_PER_NODE
  [ -z "$user_CORES_PER_NODE_CLUSTER" ] || wrt_scalar CORES_PER_NODE_CLUSTER
  [ -z "$user_TIME_FORMAT" ] || wrt_scalar TIME_FORMAT
  [ -z "$user_RELPATHNAMES" ] || wrt_scalar RELPATHNAMES
  for var in $cluster_variables_user ; do
   eval "desc=\"\$user_USER_DESCRIPTION_$var\" ;\
    allow=\"\$user_USER_ALLOWED_$var\" ;\
    deflt=\"\$user_USER_DEFAULT_$var\""
   echo "#-! USER.DESCRIPTION.$var: $desc" >> "$arch_file"
   [ -z "$deflt" ] || echo "#-! USER.DEFAULT.$var: $deflt" >> "$arch_file"
   [ -z "$allow" ] || echo "#-! USER.ALLOWED.$var: $allow" >> "$arch_file"
  done
  case "$limits_type" in
  ncore)
   wrt_scalar MAX_NCORE ;;
  indep)
   wrt_scalar MAX_NCORE
   wrt_scalar MAX_WALLTIME ;;
  queue_ncore)
   cat >> "$arch_file" <<___EOF
#-! *MAX_NCORE:
#-!  case "&USER.QUEUE&" in
___EOF
   set -- $user_MAX_NCORE_per_queue
   for queue in $user_USER_ALLOWED_QUEUE ; do
    cat >> "$arch_file" <<____EOF
#-!  $queue)
#-!   echo "$1" ;;
____EOF
    shift
   done
   wrt "#-!  esac" ;;
  queue)
   cat >> "$arch_file" <<___EOF
#-! *MAX_NCORE:
#-!  case "&USER.QUEUE&" in
___EOF
   set -- $user_MAX_NCORE_per_queue
   for queue in $user_USER_ALLOWED_QUEUE ; do
    cat >> "$arch_file" <<____EOF
#-!  $queue)
#-!   echo "$1" ;;
____EOF
    shift
   done
   wrt "#-!  esac"
   cat >> "$arch_file" <<___EOF
#-! *MAX_WALLTIME:
#-!  case "&USER.QUEUE&" in
___EOF
   set -- $user_MAX_WALLTIME_per_queue
   for queue in $user_USER_ALLOWED_QUEUE ; do
    cat >> "$arch_file" <<____EOF
#-!  $queue)
#-!   echo "$1" ;;
____EOF
    shift
   done
   wrt "#-!  esac" ;;
  interdep)
   set -- $user_MAX_WALLTIME_interdep
   i=0 ; for ncore in $user_MAX_NCORE_interdep ; do i=$((i+1))
    if ((i==1)) ; then
    cat >> "$arch_file" <<_____EOF
#-! *MAX_WALLTIME:
#-!  if ((&NCORE<=$ncore)) ; then
#-!   echo "$1"
_____EOF
    else
    cat >> "$arch_file" <<_____EOF
#-!  elif ((&NCORE<=$ncore)) ; then
#-!   echo "$1"
_____EOF
    fi
    if (($#==1)) ; then
     cat >> "$arch_file" <<______EOF
#-!  fi
#-! MAX_NCORE: $ncore
______EOF
    fi
    shift
   done ;;
  esac
  if ((review_queue>0)) ; then
   [ ! -z "$user_RUN_PARALLEL" ]\
    && [ "$user_RUN_PARALLEL" != "mpirun -np &NPROC& &BINARY&" ]\
    && wrt_scalar RUN_PARALLEL
   if [ ! -z "$user_SCRIPT_HEAD" ] ; then
    wrt "#-! SCRIPT_HEAD:"
    i=0 ; while ((i<${user_SCRIPT_HEAD[0]})) ; do i=$((i+1))
     wrt "#-!  ${user_SCRIPT_HEAD[$i]}"
    done
   fi
   if [ ! -z "$user_SCRIPT_RUN" ] ; then
    wrt "#-! SCRIPT_RUN:"
    i=0 ; while ((i<${user_SCRIPT_RUN[0]})) ; do i=$((i+1))
     wrt "#-!  ${user_SCRIPT_RUN[$i]}"
    done
   fi
   [ -z "$user_SUBMIT_SCRIPT" ] || wrt_scalar SUBMIT_SCRIPT
  fi ;;
 esac
 cmt="Includes" ; cmt1=1
 [ -z "$chosen_os" ] || wrt "include \$(INCBASE)/os/$chosen_os.arch"
 [ -z "$chosen_f90" ] || wrt "include \$(INCBASE)/f90/$chosen_f90.arch"
 [ -z "$chosen_cc" ] || wrt "include \$(INCBASE)/cc/$chosen_cc.arch"
 [ -z "$chosen_cxx" ] || wrt "include \$(INCBASE)/cxx/$chosen_cxx.arch"
 [ -z "$chosen_queue" ] || wrt "include \$(INCBASE)/queue/$chosen_queue.arch"
 cmt="Type definition" ; cmt1=1
 wrt "TYPE = $user_TYPE"
 cmt="Compiler names" ; cmt1=1
 clear_tags
 [ -z "$chosen_f90" ] || load_tags "f90/$chosen_f90.arch"
 [ -z "$chosen_cc" ] || load_tags "cc/$chosen_cc.arch"
 [ -z "$chosen_cxx" ] || load_tags "cxx/$chosen_cxx.arch"
 wrt_make F90
 wrt_make CC
 wrt_make CXX
 cmt="Environment" ; cmt1=1
 [ ! -z "$user_ENVIRONMENT_COMMAND" ]\
  && [ "$user_ENVIRONMENT_COMMAND" != ":" ]\
  && wrt_make ENVIRONMENT_COMMAND
 cmt="MPI library" ; cmt1=1
 wrt_make LIB_PATH
 wrt_make INCLUDE_DIR
 wrt_make LDLIBS_all
 wrt_make MPI_VERSION 2
 cmt="Automatic C flags" ; cmt1=1
 wrt_make CFLAGS_F90_INTERFACE
 wrt_make NEED_ETIME no
 wrt_make CFLAGS_ETIME
 cmt="Numerical libraries" ; cmt1=1
 wrt_make HAVE_BLAS no
 wrt_make HAVE_LAPACK no
 cmt="SHM support" ; cmt1=1
 wrt_make SUPPORT_SHM no
 [ "$user_SUPPORT_SHM" = yes ] && wrt_make CFLAGS_SHM
 if ((review_compiler>0)) ; then
  cmt="Compiler flags" ; cmt1=1
  wrt_make FFLAGS_opt
  wrt_make FFLAGS_debug
  wrt_make FFLAGS0_libs -O0
  cmt="Features" ; cmt1=1
  wrt_make SUPPORT_OPENMP no
  if [ "$user_SUPPORT_OPENMP" = yes ] ; then
   wrt_make FFLAGS_OPENMP_yes
   wrt_make CFLAGS_OPENMP_yes
  fi
 fi
 if [ "$user_UTILS_MODE" = native ] ; then
  cmt="Native compilation" ; cmt1=1
  wrt_make UTILS_MODE
  [ ! -z "$user_ENVIRONMENT_COMMAND_NATIVE" ]\
   && [ "$user_ENVIRONMENT_COMMAND_NATIVE" != ":" ]\
   && wrt_make ENVIRONMENT_COMMAND_NATIVE
  clear_tags
  [ -z "$chosen_f90_native" ] || load_tags "f90/$chosen_f90_native.arch"
  [ -z "$chosen_cc_native" ] || load_tags "cc/$chosen_cc_native.arch"
  wrt_make_native F90
  wrt_make_native CC
  wrt_make_native FFLAGS_opt
  wrt_make_native FFLAGS0_libs -O0
  wrt_make_native CFLAGS_F90_INTERFACE
  wrt_make_native HAVE_BLAS no
  wrt_make_native HAVE_LAPACK no
 fi
 # Additional data to help integrate into CASINO:
 cmt="Integration info (can be safely deleted)" ; cmt1=1
 wrt "### KERNEL: $host_KERNEL"
 wrt "### OS: $host_OS"
 wrt "### ARCH: $host_ARCH"
 wrt "### DISTRIBUTION: $host_DISTRIBUTION"
 wrt "### HOSTNAME: $host_HOSTNAME"
 wrt "### BASH_VERSION: $BASH_VERSION"
 for c in F90 CC CXX ; do
  eval "comp=\"\$user_$c\" ; chosen=\"\$chosen_$(uncap $c)\""
  [ -z "$comp" ] && continue
  [ -z "$chosen" ] || continue
  for flag in -version --version -V -v -h -help --help ; do
   wrt "### $c: output of $comp $flag:"
   {
    while : ; do
     IFS=""
     read -r line || break
     IFS="$IFS_save"
     set -- $line
     if [ "$1" = ArchInfoExitStatus ] ; then
      wrt "### $c: exit status of $comp $flag: $2"
     else
      wrt "###### $line"
     fi
    done
    IFS="$IFS_save"
   } < <(
    bash < <(
     echo "tmp_dir=\$(mkdir -d 2> /dev/null)"
     echo "cd \"\$tmp_dir\" >& /dev/null"
     echo "$user_ENVIRONMENT_COMMAND >& /dev/null"
     echo "$comp $flag </dev/null 2>&1"
     echo "echo \"ArchInfoExitStatus $?\""
     echo "cd >& /dev/null"
     echo "rm -rf \"\$tmp_dir\" >& /dev/null"
    ) 2> /dev/null
    (exit 0) # bypass bash fd leak (v3.2 - v4.1)
   )
   (exit 0) # bypass bash fd leak (v3.2 - v4.1)
  done
 done
 if [ -e "/proc/cpuinfo" ] ; then
  {
   wrt "### ARCH: first block of /proc/cpuinfo:"
   while : ; do
    IFS=""
    read -r line || break
    IFS="$IFS_save"
    wrt "###### $line"
   done
   IFS="$IFS_save"
  } < <(sed -n "0,/^$/p" /proc/cpuinfo)
  (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 fi
 echo -n "$el$cr"
 pretty_print 1 1 "CASINO_ARCH=$user_CASINO_ARCH defined."
 echo
 if [ ! -z "$cluster_variables_internal" ] ; then
  pretty_print 1 1 "Remember to edit CASINO/arch/data/$user_CASINO_ARCH.arch\
   to define the &INTERNAL.*& variables you have defined."
  echo
 fi
 gen_created_arch=$user_CASINO_ARCH
 return 0
}

run_compiler_test() {
 # Run compiler test dumped by function '$1', echoing:
 # 1 if compilation fails
 # 3 if the binary does not run
 # 2 if output does not match argument 'output'
 # 0 if the test succeeds
 local arg_test arg_f90 arg_cc arg_cxx arg_lflags arg_fflags arg_cflags
 local arg_cxxflags arg_output arg_env arg_mpirun arg_use_mpirun
 local tmp_dir stat run_binary
 arg_test="$1" ; shift
 while (($#>0)) ; do
  case "$1" in
  *=*) eval "arg_${1%%=*}=\"\${1#*=}\"" ;;
  *) return 1 ;;
  esac
  shift
 done
 [ -z "$arg_test" ] && echo "-1" && return
 [ -z "$arg_f90" ] && [ -z "$arg_cc" ] && [ -z "$arg_cxx" ] && echo "-1"\
  && return
 [ -z "$arg_env" ] && arg_env=":"
 [ -z "$arg_use_mpirun" ] && arg_use_mpirun=0
 tmp_dir="$(cd "$($mktemp_d_command)" >& /dev/null ; pwd)"
 [ ! -z "$arg_f90" ] && eval "dump_f90_test_$arg_test"\
  > "$tmp_dir/f90_test.f90"
 [ ! -z "$arg_cc" ] && eval "dump_cc_test_$arg_test"\
  > "$tmp_dir/cc_test.c"
 [ ! -z "$arg_cxx" ] && eval "dump_cxx_test_$arg_test"\
  > "$tmp_dir/cxx_test.cc"
 {
  echo "cd \"$tmp_dir\""
  echo "echo \"-- Compiler test: $arg_test\""
  echo "echo \"$ $arg_env\""
  echo "$arg_env"
  case "$arg_f90.$arg_cc.$arg_cxx" in
  *..)
   echo "echo \"$ $arg_f90 $arg_fflags -o test f90_test.f90 $arg_lflags\""
   echo "$arg_f90 $arg_fflags -o test f90_test.f90 $arg_lflags || exit 1" ;;
  .*.)
   echo "echo \"$ $arg_cc $arg_cflags -o test cc_test.c $arg_lflags\""
   echo "$arg_cc $arg_cflags -o test cc_test.c $arg_lflags || exit 1" ;;
  ..*)
   echo "echo \"$ $arg_cxx $arg_cxxflags -o test cxx_test.cc $arg_lflags\""
   echo "$arg_cxx $arg_cxxflags -o test cxx_test.cc $arg_lflags || exit 1" ;;
  *.*.)
   echo "echo \"$ $arg_cc $arg_cflags -c cc_test.c\""
   echo "$arg_cc $arg_cflags -c cc_test.c || exit 1"
   echo "echo \"$ $arg_f90 $arg_fflags -o test f90_test.f90 cc_test.o\
 $arg_lflags\""
   echo "$arg_f90 $arg_fflags -o test f90_test.f90 cc_test.o $arg_lflags\
    || exit 1" ;;
  *)
   echo "echo \"Can't run tests with provided set of compilers\
 $arg_f90.$arg_cc.$arg_cxx\""
   echo "exit 127" ;;
  esac
  echo "[ -x test ] || exit 1"
  if ((arg_use_mpirun==1)) && [ ! -z "$arg_mpirun" ] ; then
   var_BINARY=./test ; var_NPROC=1 ; var_NNODE=1 ; var_NCORE=1
   run_binary=$(substitute_variables "$arg_mpirun")
   unset var_BINARY var_NPROC var_NNODE var_NCORE
  else
   run_binary=./test
  fi
  if [ ! -z "$arg_output" ] ; then
   echo "echo \"$ $run_binary\""
   echo "if $run_binary ; then"
   echo " output=\$($run_binary 2>/dev/null)"
   echo "else"
   echo " exit 3"
   echo "fi"
   echo "if [ \"\$output\" = \"$arg_output\" ] ; then"
   echo " exit 0"
   echo "else"
   echo " exit 2"
   echo "fi"
  else
   echo "echo \"$ ./test\""
   echo "if ./test ; then"
   echo " exit 0"
   echo "else"
   echo " exit 3"
   echo "fi"
  fi
 } | bash >& $debug_compiler_output
 stat=$?
 rm -rf "$tmp_dir"
 return $stat
}

dump_f90_test_compile() {
 cat <<_EOF
 PROGRAM f90_test
 IMPLICIT NONE
 REAL(kind(1.0)) t1,tv(2)
 tv(1:2)=(/7.,2./)
 t1=maxval(tv)
 write(6,'(i1)')1
 END PROGRAM f90_test
_EOF
}

dump_cc_test_compile() {
 cat <<_EOF
 #include <stdio.h>
 int main(void) {
  printf("1\n");
  return 0;
 }
_EOF
}

dump_cxx_test_compile() {
 cat <<_EOF
 #include <iostream>
 int main() {
  std::cout << "1" << std::endl;
  return 0;
 }
_EOF
}

dump_f90_test_internal_etime() {
 cat <<_EOF
 PROGRAM f90_test
 IMPLICIT NONE
 REAL(kind(1.0)) t1,tv(2)
 INTERFACE
  REAL(kind(1.0)) FUNCTION etime(user_and_system_time)
   REAL(kind(1.0)),INTENT(inout) :: user_and_system_time(2)
  END FUNCTION etime
 END INTERFACE
 t1=etime(tv)
 write(6,'(i1)')1
 END PROGRAM f90_test
_EOF
}

dump_f90_test_blas_lapack() {
 cat <<_EOF
 PROGRAM f90_test
 INTEGER,PARAMETER :: dp=kind(1.d0)
 REAL(dp) a(3,3),b(3,3)
 INTEGER piv(3),ierr
 INTERFACE
  SUBROUTINE dcopy(n,dx,incx,dy,incy)
   IMPLICIT NONE
   INTEGER,INTENT(in) :: incx,incy,n
   REAL(kind(1.d0)),INTENT(in) :: dx(*)
   REAL(kind(1.d0)),INTENT(out) :: dy(*)
  END SUBROUTINE dcopy
  SUBROUTINE dgetrf(m,n,a,lda,ipiv,info)
   IMPLICIT NONE
   INTEGER,INTENT(in) :: m,n,lda
   REAL(kind(1.d0)),INTENT(inout) :: a(lda,*)
   INTEGER,INTENT(out) :: info,ipiv(*)
  END SUBROUTINE dgetrf
 END INTERFACE
 a(1:3,1)=(/2.d0,0.d0,0.d0/)
 a(1:3,2)=(/0.d0,2.d0,0.d0/)
 a(1:3,3)=(/0.d0,0.d0,2.d0/)
 call dcopy(9,a(1,1),1,b(1,1),1)
 call dgetrf(3,3,b,3,piv,ierr)
 call dgetrs('N',3,3,b,3,piv,a,3,ierr)
 write(6,'(i1)')nint(sum(a)/3.d0)
 END PROGRAM
_EOF
}

dump_f90_test_mpi() {
 cat <<_EOF
 PROGRAM f90_test
 IMPLICIT NONE
 INCLUDE 'mpif.h'
 INTEGER t1,nnodes,my_node,ierror
 t1=mpi_version
 call mpi_init(ierror)
 call mpi_comm_size(mpi_comm_world,nnodes,ierror)
 call mpi_comm_rank(mpi_comm_world,my_node,ierror)
 if(my_node==nnodes-1)write(6,'(i1)')mpi_version
 call mpi_finalize(ierror)
 END PROGRAM f90_test
_EOF
}

dump_f90_test_C_interface() {
 cat <<_EOF
 PROGRAM f90_test
 IMPLICIT NONE
 INTEGER a,b
 INTERFACE
  SUBROUTINE ftest(a,b)
  INTEGER a,b
  END SUBROUTINE ftest
 END INTERFACE
 a=1
 call ftest(a,b)
 write(6,'(i1)')b
 END PROGRAM f90_test
_EOF
}

dump_cc_test_C_interface() {
 cat <<_EOF
 #ifndef FTEST
   #ifdef F90_DOUBLE_UNDERSCORE
     #ifdef F90_CAPITALS
       #define FTEST FTEST__
     #else
       #define FTEST ftest__
     #endif
   #else
     #ifdef F90_NO_UNDERSCORE
       #ifdef F90_CAPITALS
         #define FTEST FTEST
       #else
         #define FTEST ftest
       #endif
     #else
       #ifdef F90_CAPITALS
         #define FTEST FTEST_
       #else
         #define FTEST ftest_
       #endif
     #endif
   #endif
 #endif

 void FTEST(int *x, int *y) {
  *y = *x;
 }
_EOF
}

dump_f90_test_external_etime() {
 cat <<_EOF
 PROGRAM f90_test
 IMPLICIT NONE
 REAL(kind(1.0)) t1,tv(2)
 INTERFACE
  REAL(kind(1.0)) FUNCTION etime(user_and_system_time)
   REAL(kind(1.0)),INTENT(inout) :: user_and_system_time(2)
  END FUNCTION etime
 END INTERFACE
 t1=etime(tv)
 write(6,'(i1)')1
 END PROGRAM f90_test
_EOF
}

dump_cc_test_external_etime() {
 cat <<_EOF
 #ifndef ETIME
   #ifdef F90_DOUBLE_UNDERSCORE
     #ifdef F90_CAPITALS
       #define ETIME ETIME__
       #define DTIME DTIME__
     #else
       #define ETIME etime__
       #define DTIME dtime__
     #endif
   #else
     #ifdef F90_NO_UNDERSCORE
       #ifdef F90_CAPITALS
         #define ETIME ETIME
         #define DTIME DTIME
       #else
         #define ETIME etime
         #define DTIME dtime
       #endif
     #else
       #ifdef F90_CAPITALS
         #define ETIME ETIME_
         #define DTIME DTIME_
       #else
         #define ETIME etime_
         #define DTIME dtime_
       #endif
     #endif
   #endif
 #endif
 #ifdef CASINO_NO_ETIME
   float ETIME(float*);
   float DTIME(float*);
   float ETIME(float *t) {
     t[0]=(float)0;
     t[1]=(float)0;
     return 0;
   }
   float DTIME(float *t) {
     t[0]=(float)0;
     t[1]=(float)0;
     return 0;
   }
 #else
   #include<sys/times.h>
   #include<time.h>
   #ifdef CASINO_T3E_MODE
     float ETIME(double [2]);
     float DTIME(double [2]);
     float ETIME(double t[2]) {
       struct tms buffer;
       times(&buffer);
       t[0]=(float)buffer.tms_utime/CLK_TCK;
       t[1]=(float)buffer.tms_stime/CLK_TCK;
       return t[0]+t[1];
     }
     float DTIME(double t[2]) {
       static clock_t old_values[2] = { 0, 0};
       struct tms buffer;
       times(&buffer);
       t[0]=(float)(buffer.tms_utime-old_values[0])/CLK_TCK;
       t[1]=(float)(buffer.tms_stime-old_values[1])/CLK_TCK;
       old_values[0]=buffer.tms_utime;
       old_values[1]=buffer.tms_stime;
       return t[0]+t[1];
     }
   #else
     #ifndef CLK_TCK
       #include <bits/types.h>
       extern long int __sysconf (int);
       #define CLK_TCK ((__clock_t) __sysconf (2))
     #endif
     float ETIME(float*);
     float DTIME(float*);
     float ETIME(float *t) {
       struct tms buffer;
       times(&buffer);
       t[0]=(float)buffer.tms_utime/CLK_TCK;
       t[1]=(float)buffer.tms_stime/CLK_TCK;
       return t[0]+t[1];
     }
     float DTIME(float *t) {
       static clock_t old_values[2] = { 0, 0};
       struct tms buffer;
       times(&buffer);
       t[0]=(float)(buffer.tms_utime-old_values[0])/CLK_TCK;
       t[1]=(float)(buffer.tms_stime-old_values[1])/CLK_TCK;
       old_values[0]=buffer.tms_utime;
       old_values[1]=buffer.tms_stime;
       return t[0]+t[1];
     }
   #endif
 #endif
_EOF
}

dump_f90_test_shm() {
 cat <<_EOF
 PROGRAM f90_test
 INCLUDE 'mpif.h'
 INTEGER,PARAMETER :: dp=kind(1.d0)
 INTEGER,PARAMETER :: i64=selected_int_kind(15)
 INTEGER,PARAMETER :: sz=1024,itst=512
 INTEGER ret,it1,i,nnodes,my_node,ierror,nsmps,nnpsmp
 INTEGER,ALLOCATABLE :: smp_masters(:),smp_nodes(:)
 INTEGER(mpi_address_kind) l_a
 REAL(dp),POINTER :: a(:)=>null()
 REAL(dp) aa(sz)
 POINTER (p_aa,aa) ! NON-STANDARD

 EXTERNAL alloc_shm
 EXTERNAL dealloc_shm
 EXTERNAL get_smp_list

 call mpi_init(ierror)
 call mpi_comm_size(mpi_comm_world,nnodes,ierror)
 call mpi_comm_rank(mpi_comm_world,my_node,ierror)

 allocate(smp_masters(1),smp_nodes(1))
 call get_smp_list(nsmps,smp_masters,nnpsmp,smp_nodes)

 call alloc_shm(p_aa,int(sz,i64),MPI_DOUBLE_PRECISION,ret)
 call cray_ptr_to_f90_ptr(aa)
 a(1)=1.d0
 do i=2,sz
  a(i)=a(i-1)**2
 enddo
 it1=nint(a(itst))
 call geta(a)
 call dealloc_shm(l_a,size(a),MPI_DOUBLE_PRECISION)
 nullify(a)

 write(6,'(i1)')it1

 call mpi_finalize(ierror)

 CONTAINS

  SUBROUTINE cray_ptr_to_f90_ptr(aa)
  IMPLICIT NONE
  REAL(dp),TARGET :: aa(sz)
  a=>aa
  END SUBROUTINE cray_ptr_to_f90_ptr

  SUBROUTINE geta(a)
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: a(:)
  l_a=loc(a) ! NON-STANDARD
  END SUBROUTINE geta

 END PROGRAM f90_test
_EOF
}

dump_cc_test_shm() {
 cat <<'_EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#ifdef SHM_SYSV
  #include <sys/ipc.h>
  #include <sys/shm.h>
#else
  #ifdef SHM_POSIX
    #include <unistd.h>
    #include <sys/mman.h>
    #include <fcntl.h>
  #endif
#endif
#include <errno.h>
#include <mpi.h>
#ifdef _OPEMMP
  #include <omp.h>
#endif
#ifndef ALLOC_SHM
  #ifdef F90_DOUBLE_UNDERSCORE
    #ifdef F90_CAPITALS
      #define ALLOC_SHM ALLOC_SHM__
      #define DEALLOC_SHM DEALLOC_SHM__
      #define GET_SMP_LIST GET_SMP_LIST__
    #else
      #define ALLOC_SHM alloc_shm__
      #define DEALLOC_SHM dealloc_shm__
      #define GET_SMP_LIST get_smp_list__
    #endif
  #else
    #ifdef F90_NO_UNDERSCORE
      #ifdef F90_CAPITALS
        #define ALLOC_SHM ALLOC_SHM
        #define DEALLOC_SHM DEALLOC_SHM
        #define GET_SMP_LIST GET_SMP_LIST
      #else
        #define ALLOC_SHM alloc_shm
        #define DEALLOC_SHM dealloc_shm
        #define GET_SMP_LIST get_smp_list
      #endif
    #else
      #ifdef F90_CAPITALS
        #define ALLOC_SHM ALLOC_SHM_
        #define DEALLOC_SHM DEALLOC_SHM_
        #define GET_SMP_LIST GET_SMP_LIST_
      #else
        #define ALLOC_SHM alloc_shm_
        #define DEALLOC_SHM dealloc_shm_
        #define GET_SMP_LIST get_smp_list_
      #endif
    #endif
  #endif
#endif
static MPI_Comm node_comm;
static MPI_Group group_world, node_group;
static int shm_debug=0;
static int mype, npes, mast, err;
static int numablk=0,nthreads=0;
void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes);
void ALLOC_SHM(void **ptr, int64_t *nelem, MPI_Fint *ftype, int *ret);
void DEALLOC_SHM(void **shm, int64_t *nelem, MPI_Fint *ftype);
#ifdef SHM_SYSV
  static void alloc_shm_sysv(void **ptr, size_t size, int *ret);
  static void dealloc_shm_sysv(void **shm);
#else
#ifdef SHM_POSIX
  static void alloc_shm_posix(void **ptr, size_t size, int *ret);
  static void dealloc_shm_posix(void **shm, int64_t *nelem, MPI_Datatype type);
#endif
#endif
static void set_shm_numablock();
#define NUMATAG_WIDTH 5
#ifdef SHM_SYSV
  void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
    int nlen_processor, nlen;
    int i, n, found, nn=0, np=0;
    char nodnam[MPI_MAX_PROCESSOR_NAME+NUMATAG_WIDTH], *cnidlst, *tc1, *tc2;
    char string[20];
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Get_processor_name(nodnam, &nlen_processor);
    MPI_Allreduce(&nlen_processor, &nlen, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    set_shm_numablock();
    if ( numablk > 0 ) {
      if ( nthreads == 0 ) {
        n = mype/numablk;
      } else {
        n = mype/(numablk/nthreads);
      }
      nlen += NUMATAG_WIDTH;
      sprintf(string,"%-5i",n);
      strncat(nodnam,string,NUMATAG_WIDTH);
    }
    cnidlst = (char *) malloc((npes*nlen+1)*sizeof(char));
    MPI_Allgather(&nodnam, nlen, MPI_CHAR, cnidlst, nlen, MPI_CHAR,
       MPI_COMM_WORLD);
    nn = np = 0;
    for (i=0; i<npes; i++){
      tc1=cnidlst+i*nlen;
      if( strncmp(nodnam,tc1,nlen) == 0 ) nodpes[np++]=i;
      found = 0;
      for (n=0; n<nn; n++) {
        tc2=cnidlst+nodlst[n]*nlen;
        if ( strncmp(tc1,tc2,nlen) == 0 )  found = 1;
      }
      if (!found) nodlst[nn++] = i;
    }
    *nnod = nn;
    *nppn = np;
    mast = nodpes[0];
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);
    MPI_Group_incl(group_world,np,nodpes,&node_group);
    MPI_Comm_create(MPI_COMM_WORLD,node_group,&node_comm);
    free(cnidlst);
  }
#else
  #ifdef SHM_POSIX
    static char posix_shm_fn[20]="/casino.shm";
    void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
      int i, j, is, ie, n, fd, shift, found, nn=0, np=0, tag0=-1;
      char string[20];
      const char fname[]="/casino_shm_tmp";
      off_t size;
      int *nd_shm, *tag;
      MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      MPI_Comm_size(MPI_COMM_WORLD, &npes);
      tag = (int *)malloc(npes*sizeof(int));
      for(i=0;i<npes;i++) tag[i]=-1;
      size=(off_t)npes*(off_t)sizeof(int);
      sprintf(string,"pe %d : shm_open",mype);
      if ((fd=shm_open(fname,(O_CREAT|O_RDWR), 0666)) < 0) perror(string);
      sprintf(string,"pe %d : ftruncate",mype);
      if (ftruncate(fd,size) < 0) perror(string);
      sprintf(string,"pe %d : mmap",mype);
      if ((nd_shm = (int *)mmap(NULL, size, (PROT_READ|PROT_WRITE),
         MAP_SHARED, fd, 0))< 0) perror(string);
      MPI_Barrier(MPI_COMM_WORLD);
      for(i=0;i<npes;i++) nd_shm[i]=-1;
      MPI_Barrier(MPI_COMM_WORLD);
      nd_shm[mype]=0;
      set_shm_numablock();
      MPI_Barrier(MPI_COMM_WORLD);
      if(numablk == 0) {
        for(i=0; i<npes; i++) {
          if(nd_shm[i] == 0) { tag0=i; break; }
        }
        if (mype == tag0) {
          sprintf(string,"pe %d : shm_unlink",mype);
          if (shm_unlink(fname) < 0) perror(string);
        }
      } else {
        nn=0; found=1;
        for(i=0; i<npes; i++) {
          if(nd_shm[i] == 0) {
            nn++;
            if(found) { found=0; is=i; }
            ie=i;
          }
          if (mype == is) {
            sprintf(string,"pe %d : shm_unlink",mype);
            if (shm_unlink(fname) < 0) perror(string);
          }
        }
        if( (ie-is+1) != nn ) {
          fprintf(stderr,
             "pe %d : ranks are not sequentially placed on the NUMA node.\n"
             "In nd_shm we got start %d end %d total %d \n  "
             "Try to use MPI environment variable place  sequential ranks\n"
             "on the same node. Exiting!", mype, is, ie, nn);
          MPI_Barrier(MPI_COMM_WORLD);
          MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
          shift=numablk;
          if (nthreads>0) shift=numablk/nthreads;
          tag0 = is+(mype-is)/shift;
        }
      }
      MPI_Allgather(&tag0,1,MPI_INT,tag,1,MPI_INT,MPI_COMM_WORLD);
      nn=np=0;
      for(i=0;i<npes;i++){
        if(tag[mype] == tag[i]) nodpes[np++]=i;
        found=0;
        for(n=0;n<nn;n++) {
          if(tag[i] == tag[nodlst[n]]) found = 1;
        }
        if (!found) nodlst[nn++]=i;
      }
      *nnod = nn;
      *nppn = np;
      mast = nodpes[0];
      MPI_Comm_group(MPI_COMM_WORLD, &group_world);
      MPI_Group_incl(group_world, np, nodpes, &node_group);
      MPI_Comm_create(MPI_COMM_WORLD, node_group, &node_comm);
      if (numablk>0) {
        sprintf(string,"%-5i",(mype/shift));
        strncat(posix_shm_fn,string,NUMATAG_WIDTH);
      }
      sprintf(string,"pe %d : munmap size %ld",mype,size);
      if (munmap((void *)nd_shm, size) < 0) perror(string);
      free(tag);
    }
  #else
    void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
      int mype;
      MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      if (mype==0) printf(
         "It looks that the Makefile flag SHM_VERSION_ contains\n",
         "neither SHM_SYS or SHM_POSIX.  Check your Makefile flags.\n",
         "Quitting ...\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,-1);
    }
  #endif
#endif
void ALLOC_SHM(void **ptr, int64_t *nelem, MPI_Fint *ftype, int *ret) {
  int typsiz;
  size_t size;
  MPI_Datatype type;
  err = 0;
  type = MPI_Type_f2c( *ftype);
  MPI_Type_size(type, &typsiz);
  size = (size_t)*nelem * (size_t)typsiz;
  #ifdef SHM_SYSV
    alloc_shm_sysv(ptr, size, ret);
  #endif
  #ifdef SHM_POSIX
    alloc_shm_posix(ptr, size, ret);
  #endif
  if (*ret) *ret = err;
}
void DEALLOC_SHM(void **shm, int64_t *nelem, MPI_Fint *ftype) {
  MPI_Datatype type;
  type = MPI_Type_f2c( *ftype);
  #ifdef SHM_SYSV
    dealloc_shm_sysv(shm);
  #else
    #ifdef SHM_POSIX
      dealloc_shm_posix(shm, nelem, type);
    #endif
  #endif
}
#ifdef SHM_SYSV
  static void alloc_shm_sysv(void **ptr, size_t size, int *ret) {
    char string[20];
    int shmid;
    void *shm;
    struct shmid_ds ds;
    if (mype == mast) {
      if ((shmid = shmget(IPC_PRIVATE, size, 0666)) < 0) {
        sprintf(string,"pe %d: shmget",mype);
        perror(string);
        if(errno==EINVAL) printf(
           "EINVAL: requested size of shared memory segment too large (%lx).\n"
           "Ask administrator to increase system-wide limit.\n"
           "E.g.: on Linux use 'sysctl kernel.shmmax=$((1<<63))'\n", size);
        else if(errno==ENOMEM) printf(
           "ENOMEM: memory allocation for control structures failed.\n");
        else if(errno==ENOSPC) printf(
           "ENOSPC: either the system-wide limit for shared memory (SHMALL)\n"
           "or the maximum number of SHM ids (SHMMNI) has been reached.\n"
           "Maybe some other process is not cleaning up properly?\n"
           "Clean up, reboot or increase limits (e.g., on Linux use sysctl)\n");
        else printf(
           "Have hit an unidentified error when calling shmget.");
        if (*ret) err++; else exit(1);
      }
      MPI_Bcast(&shmid, 1, MPI_INT, 0, node_comm);
      shm = shmat(shmid, NULL, 0);
      if (shm == (void *) -1) {
        sprintf(string,"pe %d: shmat",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }
      if (shmctl(shmid, IPC_RMID, &ds) < 0) {
        sprintf(string,"pe %d: shmget",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }
    } else {
      MPI_Bcast(&shmid, 1, MPI_INT, 0, node_comm);
      shm = shmat(shmid, NULL, 0);
      if (shm == (void *) -1) {
        sprintf(string,"pe %d: shmat",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }
    }
    MPI_Barrier(node_comm);
    *ptr = shm;
  }
  static void dealloc_shm_sysv(void **ptr) {
    char string[20];
    if (shmdt(*ptr) < 0) {
      sprintf(string,"pe %d: shmdt",mype);
      perror(string);
      exit(1);
    }
  }
#else
  #ifdef SHM_POSIX
    static void alloc_shm_posix(void **ptr, size_t size, int *ret) {
      void *shm=NULL;
      char string[20];
      int fd,node_id;
      sprintf(string,"pe %d: shm_open",mype);
      if ((fd = shm_open(posix_shm_fn, (O_CREAT|O_RDWR), 0666)) < 0) {
        perror(string);
        if (*ret) err++; else exit(1);
      }
      sprintf(string,"pe %d: ftruncate",mype);
      if (ftruncate(fd, (off_t)size) < 0) {
        perror(string);
        if (*ret) err++; else exit(1);
      }
      sprintf(string,"pe %d: mmap",mype);
      if ((shm = mmap(NULL, size, (PROT_READ|PROT_WRITE), MAP_SHARED, fd, 0))
         == (void *) -1) {
        perror(string);
        if (errno==ENOMEM) printf(
           "ENOMEM: requested shared memory size exceeds system limit\n");
        else if (errno==EBADF) printf(
           "EBADF: wrong value for fd argument. \n");
        else if (errno==EACCES) printf(
           "EACCES: wrong access flags.\n");
        else printf(
           "Have hit an unidentified error when calling mmap.");
        if (*ret) err++; else exit(1);
      }
      MPI_Barrier(node_comm);
      MPI_Comm_rank(node_comm, &node_id);
      if (node_id == 0) {
        sprintf(string,"pe %d: shm_unlink",mype);
        if (shm_unlink(posix_shm_fn) < 0 ){
          perror(string);
          if (*ret) err++; else exit(1);
        }
      }
      *ptr = shm;
    }
    static void dealloc_shm_posix(void **shm, int64_t *nelem,
       MPI_Datatype type) {
      char string[20];
      int typsiz;
      size_t size;
      MPI_Type_size(type, &typsiz);
      size = (size_t)*nelem * (size_t)typsiz;
      sprintf(string,"pe %d : munmap size %ld", mype, size);
      if (munmap(*shm, size)<0) {
        perror(string); exit(1);
      }
    }
  #endif
#endif
void set_shm_numablock() {
  char *nb;
  nb = getenv("CASINO_NUMABLK");
  if ( nb != NULL ) numablk=atoi(nb);
  #ifdef _OPENMP
    #pragma omp parallel
    {
      #pragma omp master
      nthreads=omp_get_num_threads();
    }
  #endif
}
_EOF
}

abort() {
 # Make sure we leave no background processes
 echo "[ABORT] Killing all subprocesses...$el"
 [ -e "$tmpdir/.tmp_$script_name" ] && rm -f "$tmpdir/.tmp_$script_name"
 killall -9 $script_name
 exit 1
}
############################ END FUNCTIONS ############################

return >& /dev/null || :

init_tput
[ -e "$(dirname $0)/taglib.sh" ] || errstop "$(dirname $0)/taglib.sh library\
 not found."
source "$(dirname $0)/taglib.sh"
initialize
parse_cmd_line "$@"
post_config
# Go into directory where .arch files live.
[ -d "$(dirname $0)/data" ] || errstop "Directory $(dirname $0)/data not found"
cd "$(dirname $0)/data"
build_file_list
# Run selected action
case "$run_mode" in
auto)
 get_host_params
 trap abort INT QUIT TERM
 autodetect_arch ;;
list) list_arch $span ;;
gen)
 get_host_params
 gen_arch ;;
esac
