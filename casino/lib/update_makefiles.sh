#!/bin/bash
#----------------------------------#
# Script to re-generate the LAPACK #
# and BLAS Makefiles in CASINO/lib #
# PLR 12.2008                      #
#----------------------------------#
set +u
export LANG="POSIX"

# Tunable options
# ===============
# 1. Directories. The ones below are for TCM.
lapack_dir=~pl275/software/lapack/src/lapack
blas_dir=~pl275/software/lapack/src/blas
# 2. Files to be compiled without optimization. In the list below, some files
#  appear because they are not optimized in the original LAPACK distribution,
#  while others give problems with certain compilers in CASINO. Be careful
#  if removing any of these from the list.
noopt_lapack="slamch.f dlamch.f slaruv.f dlaruv.f"
noopt_blas=""

################################# FUNCTIONS #################################
function errstop {
 [ ! -z "$lapack_mktemp" ] && [ -e "$lapack_mktemp" ] && rm -f "$lapack_mktemp"
 [ ! -z "$blas_mktemp" ] && [ -e "$blas_mktemp" ] && rm -f "$blas_mktemp"
 echo ; echo "$*" ; echo ; exit 1
}

function nfield { echo $# ; }

function cap {
 # Turn lower case into upper case in $1.
 local i string out="" n c
 i=0 ; string="$1" ; n=${#string}
 while ((i<n)) ; do c="${string:$i:1}" ; i=$((i+1))
  case "$c" in
  a) c=A ;; b) c=B ;; c) c=C ;; d) c=D ;; e) c=E ;; f) c=F ;; g) c=G ;;
  h) c=H ;; i) c=I ;; j) c=J ;; k) c=K ;; l) c=L ;; m) c=M ;; n) c=N ;;
  o) c=O ;; p) c=P ;; q) c=Q ;; r) c=R ;; s) c=S ;; t) c=T ;; u) c=U ;;
  v) c=V ;; w) c=W ;; x) c=X ;; y) c=Y ;; z) c=Z ;;
  esac
  out="$out$c"
 done
 echo "$out"
}

function uncap {
 # Turn upper case into lower case in $1.
 local i string out="" n c
 i=0 ; string="$1" ; n=${#string}
 while ((i<n)) ; do c="${string:$i:1}" ; i=$((i+1))
  case "$c" in
  A) c=a ;; B) c=b ;; C) c=c ;; D) c=d ;; E) c=e ;; F) c=f ;; G) c=g ;;
  H) c=h ;; I) c=i ;; J) c=j ;; K) c=k ;; L) c=l ;; M) c=m ;; N) c=n ;;
  O) c=o ;; P) c=p ;; Q) c=q ;; R) c=r ;; S) c=s ;; T) c=t ;; U) c=u ;;
  V) c=v ;; W) c=w ;; X) c=x ;; Y) c=y ;; Z) c=z ;;
  esac
  out="$out$c"
 done
 echo "$out"
}

function rem_list {
 # Remove item $1 from $2-$n
 local item out=""
 for item in ${@:2} ; do [ "$item" = "$1" ] || out="$out $item" ; done
 echo "$out"
}

function in_line {
 # Like which_field, but only set return value, no other output
 local str="$1"
 while (($#>1)) ; do shift ; [ "$str" = "$1" ] && return 0 ; done
 return 1
}

function check_routine_in_file {
 # Check whether routine/function $1 is in file $2. Return status indicates
 # success.
 local r="$1" df="$2"
 grep -qE "^ *SUBROUTINE *$r|^ *FUNCTION *$r|^ *REAL(\*[0-9]*)? *FUNCTION *$r|\
^ *COMPLEX(\*[0-9]*)? *FUNCTION *$r|^ *DOUBLE *PRECISION *FUNCTION *$r|\
^ *DOUBLE *COMPLEX *FUNCTION *$r|^ *INTEGER *FUNCTION *$r|\
^ *LOGICAL *FUNCTION *$r|^ *CHARACTER *FUNCTION *$r" "$df" >& /dev/null
}

function locate_routine {
 # Locate a routine $1 in $2/*.f , $3/*.f ... $n/*.f , and return
 # directory/file.f on stdout, or nothing if not found.
 local r i df f
 r="$1" ; shift
 # Do a quick search
 r=$(cap $r) ; f=$(uncap $r).f
 i=0 ; while ((i<$#)) ;  do i=$((i+1))
  df="${@:$i:1}/$f"
  [ -e "$df" ] && check_routine_in_file "$r" "$df" && { echo "$df" ; return ; }
 done
 # Do a full search
 i=0 ; while ((i<$#)) ; do i=$((i+1))
  d="${@:$i:1}"
  { while : ; do
   read df || break
   f="${df##*/}" ; df="$d/$f"
   [ -e "$df" ] && check_routine_in_file "$r" "$df" && { echo "$df" ; return ; }
  done ; } < <(find "$d" -name "*.f")
 done
}

function copy_mfile_chunk {
 # Copy chunk of the Makefile, or give an error if appropriate
 local mkin="$1" mkout="$2" mark1="$3" ; mark2="$4"
 [ ! -z "$mark1" ] && echo "$mark1" >> $mkout
 if ! cat_partial "$mark1" "$mark2" $mkin >> $mkout ; then
  [ -z "$mark1" ] && errstop "Line '$mark2' not found in $mkin."
  [ -z "$mark2" ] && errstop "Line '$mark1' not found in $mkin."
  errstop "Lines '$mark1' and/or '$mark2' not found in $mkin."
 fi
 [ ! -z "$mark2" ] && echo "$mark2" >> $mkout
}

function cat_partial {
 # Dump the contents of file $3 between lines containing "$1" and "$2".
 local string0="$1" string1="$2" file=$3 saved_IFS="$IFS" line mode=0
 export IFS=""
 { [ -z "$string0" ] && mode=1
 while ((mode==0)) ; do
  read -r line || break ; [ "$line" = "$string0" ] && mode=1
 done
 ((mode==0)) && { export IFS="$saved_IFS" ; return 1 ; }
 [ -z "$string1" ] && mode=2
 while ((mode==1)) ; do
  read -r line || break ; [ "$line" = "$string1" ] && mode=0 || echo "$line"
 done
 ((mode!=2)) && { export IFS="$saved_IFS" ; return $mode ; }
 while : ; do
  read -r line || break ; echo "$line"
 done ; } <$file
 export IFS="$saved_IFS" ; return 0
}
############################### END FUNCTIONS ###############################

# Print header
echo "LAPACK source picker"
echo "===================="
echo "Script to find what source files you need in order to use a given set"
echo "of LAPACK and BLAS routines."
echo

# Set up
[[ "$TERM" == xterm-* ]] && export TERM=xterm
el=$(tput el) ; cr=$(tput cr)
while : ; do
 [ -d "$blas_dir" ] && break
 echo "Could not find BLAS sources directory $blas_dir."
 echo
 echo "In the original LAPACK 3.1.1 tarball the BLAS sources are in BLAS/SRC/"
 echo "Enter directory below (empty to quit):"
 echo -n "> " ; read blas_dir
 echo
 [ -z "$blas_dir" ] && errstop "Quitting."
done
while : ; do
 [ -d "$lapack_dir" ] && break
 echo "Could not find LAPACK sources directory $lapack_dir."
 echo
 echo "In the original LAPACK 3.1.1 tarball the LAPACK sources are in SRC/,"
 echo "but bear in mind that there are some routines are under INSTALL/ too:"
 echo "if you get a 'routine not found' error later on you will need to browse"
 echo "the INSTALL/ directory, find the routine and copy it into SRC/."
 echo
 echo "Enter directory below (empty to quit):"
 echo -n "> " ; read lapack_dir
 echo
 [ -z "$lapack_dir" ] && errstop "Quitting."
done
blas_dir=$(cd "$blas_dir" ; pwd) ; blas_dir=${blas_dir%/}
lapack_dir=$(cd "$lapack_dir" ; pwd) ; lapack_dir=${lapack_dir%/}
[ "$blas_dir" = "$lapack_dir" ] &&\
 errstop "BLAS and LAPACK sources should not be in the same directory."
echo "BLAS sources directory: $blas_dir"
echo "LAPACK sources directory: $lapack_dir"
echo

# Set up CASINO directory
casino_blas="./BLAS"
casino_lapack="./LAPACK"
[ -d "$casino_blas" ] && [ -d "$casino_lapack" ]\
 || errstop "This script must be run under CASINO/lib ."
lapack_mk="$casino_lapack/Makefile"
blas_mk="$casino_blas/Makefile"
[ -f "$lapack_mk" ] && [ -f "$blas_mk" ]\
 || errstop "This script must be run under CASINO/lib ."
lapack_mktemp="$lapack_mk.temp_$$"
blas_mktemp="$blas_mk.temp_$$"

# Read routines
echo "Enter the names of the LAPACK and BLAS routines which you plan to use"
echo "in CASINO (case-insensitive, space-separated list; empty to quit):"
echo -n "> " ; read initial_routine_list
echo
[ -z "$initial_routine_list" ] \
 && { echo "Empty list. Quitting." ; exit ; }

# Get file names for routines, and sort them into LAPACK/BLAS
lapack_rlist_clevel="" ; lapack_rlist=""
lapack_flist_clevel="" ; lapack_flist=""
blas_rlist_clevel="" ; blas_rlist=""
blas_flist_clevel="" ; blas_flist=""
set -- $initial_routine_list
for r in $* ; do
 r=$(cap $r)
 in_line $r $lapack_rlist_clevel $blas_rlist_clevel &&\
  { echo "WARNING: duplicate LAPACK entry $r" ; continue ; }
 df=$(locate_routine "$r" "$lapack_dir" "$blas_dir")
 [ -z "$df" ] && errstop "Routine $r not found."
 d=${df%/*.f} ; f=${df##$d/}
 if [ "$d" = "$lapack_dir" ] ; then
  lapack_rlist_clevel="$lapack_rlist_clevel $r"
  in_line $f $lapack_flist_clevel\
   || lapack_flist_clevel="$lapack_flist_clevel $f"
 elif [ "$d" = "$blas_dir" ] ; then
  blas_rlist_clevel="$blas_rlist_clevel $r"
  in_line $f $blas_flist_clevel\
   || blas_flist_clevel="$blas_flist_clevel $f"
 else
  errstop "Problem locating routine $r: locator returned '$df'."
 fi
done

# Go over LAPACK dependencies
ilevel=0 ; while (($(nfield $lapack_rlist_clevel)>0)) ; do ilevel=$((ilevel+1))
 nrout=$(nfield $lapack_rlist_clevel)
 nfile=$(nfield $lapack_flist_clevel)
 echo "LAPACK cascade level $ilevel contains $nrout routines in $nfile source\
 files.$el"
 # Add routines/files in current level to list
 lapack_flist="$lapack_flist $lapack_flist_clevel"
 lapack_rlist="$lapack_rlist $lapack_rlist_clevel"
 # Loop over files in current level
 flist="" ; rlist=""
 ifile=0 ; for cf in $lapack_flist_clevel ; do ifile=$((ifile+1))
  [ -e "$lapack_dir/$cf" ] || errstop "File $cf not found in LAPACK directory."
  echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf$el$cr"
  # Parse current file for EXTERNALs
  {
   while : ; do
    read line || continue 2
    set -- ${line//,/, }
    if [ "$1" = EXTERNAL ] ; then
     shift # discard 'EXTERNAL'
     # Loop over routines mentioned in EXTERNAL statement
     while (($#>0)) ; do
      r=$1
      # Flag presence of comma, which should indicate continuation line
      [ "${r:$((${#r}-1)):1}" = , ] && comma=1 || comma=0
      r=$(cap ${r%,})
      echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r$el$cr"
      if ! in_line $r $lapack_rlist $rlist \
       && ! in_line $r $blas_rlist_clevel ; then
       # Not in routine list yet, try to locate routine
       df=$(locate_routine "$r" "$lapack_dir" "$blas_dir")
       [ -z "$df" ] && errstop "Could not find routine $r."
       d=${df%/*.f} ; f=${df##$d/}
       echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r < $f$el$cr"
       if [ "$d" = "$lapack_dir" ] ; then # LAPACK routine
        echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r < $f\
 < LAPACK$el$cr"
        rlist="$rlist $r"
        in_line $f $lapack_flist $flist || flist="$flist $f"
       elif [ "$d" = "$blas_dir" ] ; then # BLAS routine
        echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r < $f\
 < BLAS$el$cr"
        blas_rlist_clevel="$blas_rlist_clevel $r"
        in_line $f $blas_flist_clevel\
         || blas_flist_clevel="$blas_flist_clevel $f"
       else # not BLAS or LAPACK: error
        errstop "Problem locating routine $r: locator returned '$df'."
       fi # location of routine
      fi # new routine
      shift
      if (($#<1)) && ((comma==1)) ; then
       read line
       set -- ${line//,/, }
       shift # discard continuation character
      fi
     done
    fi
   done
  } < "$lapack_dir/$cf"
 done
 # Dependencies at current level become new level
 lapack_flist_clevel="$flist"
 lapack_rlist_clevel="$rlist"
done
echo -n "$el"
echo "Done processing LAPACK."
echo

# Now go over BLAS dependencies
ilevel=0 ; while (($(nfield $blas_rlist_clevel)>0)) ; do ilevel=$((ilevel+1))
 nrout=$(nfield $blas_rlist_clevel)
 nfile=$(nfield $blas_flist_clevel)
 echo "BLAS cascade level $ilevel contains $nrout routines in $nfile source\
 files.$el"
 # Add routines/files in current level to list
 blas_flist="$blas_flist $blas_flist_clevel"
 blas_rlist="$blas_rlist $blas_rlist_clevel"
 # Loop over files in current level
 flist="" ; rlist=""
 ifile=0 ; for cf in $blas_flist_clevel ; do ifile=$((ifile+1))
  [ -e "$blas_dir/$cf" ] || errstop "File $cf not found in BLAS directory."
  echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf$el$cr"
  # Parse current file for EXTERNALs
  {
   while : ; do
    read line || continue 2
    set -- $line
    if [ "$1" = EXTERNAL ] ; then
     shift # discard 'EXTERNAL'
     # Loop over routines mentioned in EXTERNAL statement
     while (($#>0)) ; do
      r=$1
      # Flag presence of comma, which should indicate continuation line
      [ "${r:$((${#r}-1)):1}" = , ] && comma=1 || comma=0
      r=$(cap ${r%,})
      echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r$el$cr"
      if ! in_line $r $blas_rlist $rlist ; then
       # Not in routine list yet, try to locate routine
       df=$(locate_routine "$r" "$blas_dir")
       echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r < $f$el$cr"
       [ -z "$df" ] && errstop "Could not find routine $r."
       d=${df%/*.f} ; f=${df##$d/}
       if [ "$d" = "$blas_dir" ] ; then # BLAS routine
        echo -n "Scanning level $ilevel... ($ifile/$nfile) $cf < $r < $f\
 < BLAS$el$cr"
        rlist="$rlist $r"
        in_line $f $blas_flist $flist || flist="$flist $f"
       else # not BLAS: error
        errstop "Problem locating routine $r: locator returned '$df'."
       fi # location of routine
      fi # new routine
      shift
      if (($#<1)) && ((comma==1)) ; then
       read line
       set -- $line
       shift # discard continuation character
      fi
     done
    fi
   done
  } < "$blas_dir/$cf"
 done
 # Dependencies at current level become new level
 blas_flist_clevel="$flist"
 blas_rlist_clevel="$rlist"
done
echo -n "$el"
echo "Done processing BLAS."
echo

# Sort file lists alphabetically
lapack_flist=$( { for f in $lapack_flist ; do echo $f ; done ; } | sort )
blas_flist=$( { for f in $blas_flist ; do echo $f ; done ; } | sort )

lapack_any_action=0 ; blas_any_action=0
lapack_mk_changed=0 ; blas_mk_changed=0

# Determine what to do with LAPACK
echo "For LAPACK:"
echo
lapack_copy="" ; lapack_replace="" ; lapack_remove=""
for f in $lapack_flist ; do
 if [ -f "$casino_lapack/$f" ] ; then
  cmp "$lapack_dir/$f" "$casino_lapack/$f" >& /dev/null\
   || lapack_replace="$lapack_replace $f"
 else
  lapack_copy="$lapack_copy $f"
 fi
done
{ while : ; do
 read df || break
 f="${df##*/}" ; df="$casino_lapack/$f"
 [ -e "$df" ] || continue
 [ -e "$lapack_dir/$f" ] && in_line $f $lapack_flist && continue
 lapack_remove="$lapack_remove $f"
done ; } < <(find "$casino_lapack" -name "*.f")

# Report file actions
[ -z "$lapack_copy" ] || echo -e "New files to copy over: $lapack_copy\n"
[ -z "$lapack_replace" ] || echo -e "Files to replace: $lapack_replace\n"
[ -z "$lapack_remove" ] || echo -e "Files to remove: $lapack_remove\n"
[ -z "$lapack_copy$lapack_replace$lapack_remove" ] || lapack_any_action=1

# Build new Makefile
touch "$lapack_mktemp"
copy_mfile_chunk "$lapack_mk" "$lapack_mktemp" "" "# Object-file list"
for add in "" _NATIVE ; do
 # Optimized
 line="OBJ_OPT$add ="
 for f in $lapack_flist ; do
  in_line $f $noopt_lapack && continue
  item=" \$(OBJDIR$add)/${f%.f}.o"
  if ((${#line}+${#item}>78)) ; then
   echo "$line \\" >> "$lapack_mktemp"
   line=""
  fi
  line="$line$item"
 done
 [ -z "$line" ] || echo "$line" >> "$lapack_mktemp"
 # Not optimized
 line="OBJ_NOOPT$add ="
 for f in $lapack_flist ; do
  in_line $f $noopt_lapack || continue
  item=" \$(OBJDIR$add)/${f%.f}.o"
  if ((${#line}+${#item}>78)) ; then
   echo "$line \\" >> "$lapack_mktemp"
   line=""
  fi
  line="$line$item"
 done
 [ -z "$line" ] || echo "$line" >> "$lapack_mktemp"
done
# Rest of Makefile
copy_mfile_chunk "$lapack_mk" "$lapack_mktemp" "# End object-file list" ""

# Report Makefile changes
lapack_mk_changed=0
if ! cmp "$lapack_mk" "$lapack_mktemp" >& /dev/null ; then
 lapack_any_action=1 ; lapack_mkfile_changed=1
 echo "Makefile changes:"
 diff "$lapack_mk" "$lapack_mktemp"
 echo
else
 rm -f "$lapack_mktemp"
fi

((lapack_any_action==0)) && echo -e "No action required.\n"

# Determine what to do with BLAS
echo "For BLAS:"
echo
blas_copy="" ; blas_replace="" ; blas_remove=""
for f in $blas_flist ; do
 if [ -f "$casino_blas/$f" ] ; then
  cmp "$blas_dir/$f" "$casino_blas/$f" >& /dev/null\
   || blas_replace="$blas_replace $f"
 else
  blas_copy="$blas_copy $f"
 fi
done
{ while : ; do
 read df || break
 f="${df##*/}" ; df="$casino_blas/$f"
 [ -e "$df" ] || continue
 [ -e "$blas_dir/$f" ] && in_line $f $blas_flist && continue
 blas_remove="$blas_remove $f"
done ; } < <(find "$casino_blas" -name "*.f")

# Report file actions
[ -z "$blas_copy" ] || echo -e "New files to copy over: $blas_copy\n"
[ -z "$blas_replace" ] || echo -e "Files to replace: $blas_replace\n"
[ -z "$blas_remove" ] || echo -e "Files to remove: $blas_remove\n"
[ -z "$blas_copy$blas_replace$blas_remove" ] || blas_any_action=1

# Build new Makefile
touch "$blas_mktemp"
copy_mfile_chunk "$blas_mk" "$blas_mktemp" "" "# Object-file list"
for add in "" _NATIVE ; do
 # Optimized
 line="OBJ_OPT$add ="
 for f in $blas_flist ; do
  in_line $f $noopt_blas && continue
  item=" \$(OBJDIR$add)/${f%.f}.o"
  if ((${#line}+${#item}>78)) ; then
   echo "$line \\" >> "$blas_mktemp"
   line=""
  fi
  line="$line$item"
 done
 [ -z "$line" ] || echo "$line" >> "$blas_mktemp"
 # Not optimized
 line="OBJ_NOOPT$add ="
 for f in $blas_flist ; do
  in_line $f $noopt_blas || continue
  item=" \$(OBJDIR$add)/${f%.f}.o"
  if ((${#line}+${#item}>78)) ; then
   echo "$line \\" >> "$blas_mktemp"
   line=""
  fi
  line="$line$item"
 done
 [ -z "$line" ] || echo "$line" >> "$blas_mktemp"
done
# Rest of Makefile
copy_mfile_chunk "$blas_mk" "$blas_mktemp" "# End object-file list" ""

# Report Makefile changes
blas_mk_changed=0
if ! cmp "$blas_mk" "$blas_mktemp" >& /dev/null ; then
 blas_any_action=1 ; blas_mkfile_changed=1
 echo "Makefile changes:"
 diff "$blas_mk" "$blas_mktemp"
 echo
else
 rm -f "$blas_mktemp"
fi

((blas_any_action==0)) && echo -e "No action required.\n"

((lapack_any_action==0)) && ((blas_any_action==0)) && exit
echo "Perform these actions [Y/n]?"
echo -n "> " ; read yn
[ -z "$yn" ] && yn=Y
case "$yn" in
y*|Y*) : ;; *) errstop "Quitting." ;;
esac

# Actions on LAPACK
if ((lapack_any_action==1)) ; then
 echo -n "Updating LAPACK...$el$cr"
 for f in $lapack_copy ; do
  echo -n "Updating LAPACK... copying $f$el$cr"
  cp "$lapack_dir/$f" "$casino_lapack/$f"
 done
 for f in $lapack_replace ; do
  echo -n "Updating LAPACK... replacing $f$el$cr"
  rm -f "$casino_lapack/$f"
  cp "$lapack_dir/$f" "$casino_lapack/$f"
 done
 for f in $lapack_remove ; do
  echo -n "Updating LAPACK... removing $f$el$cr"
  rm -f "$casino_lapack/$f"
 done
 if ((lapack_mkfile_changed==1)) ; then
  echo -n "Updating LAPACK... updating Makefile$el$cr"
  mv -f "$lapack_mktemp" "$lapack_mk"
 fi
 echo "LAPACK updated.$el"
fi

# Actions on BLAS
if ((blas_any_action==1)) ; then
 echo -n "Updating BLAS...$el$cr"
 for f in $blas_copy ; do
  echo -n "Updating BLAS... copying $f$el$cr"
  cp "$blas_dir/$f" "$casino_blas/$f"
 done
 for f in $blas_replace ; do
  echo -n "Updating BLAS... replacing $f$el$cr"
  rm -f "$casino_blas/$f"
  cp "$blas_dir/$f" "$casino_blas/$f"
 done
 for f in $blas_remove ; do
  echo -n "Updating BLAS... removing $f$el$cr"
  rm -f "$casino_blas/$f"
 done
 if ((blas_mkfile_changed==1)) ; then
  echo -n "Updating BLAS... updating Makefile$el$cr"
  mv -f "$blas_mktemp" "$blas_mk"
 fi
 echo "BLAS updated.$el"
fi
