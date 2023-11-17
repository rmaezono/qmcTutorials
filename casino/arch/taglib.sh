#!/bin/bash
# Library with functions for handling the tag system in .arch files

# Don't whinge about undefined variables
set +u
# Enable extended pattern-matching features
shopt -s extglob

# Required basic functions (use host script's if already defined):
if [ "$(type -t field)" != function ] ; then
 field() { local i=$1 ; shift ; ((i>0)) && echo "${@:$i:1}" ; }
fi
if [ "$(type -t in_line)" != function ] ; then
 in_line() {
  local str="$1"
  while (($#>1)) ; do shift ; [ "$str" = "$1" ] && return 0 ; done
  return 1
 }
fi
if [ "$(type -t rem_list)" != function ] ; then
 rem_list() {
  # Remove item $1 from $2-$n
  local item out=""
  for item in ${@:2} ; do [ "$item" = "$1" ] || out="$out $item" ; done
  echo "$out"
 }
fi
if [ "$(type -t check_number_N)" != function ] ; then
 check_number_N() { [[ "$1" == +([0-9]) ]] ; }
fi
if [ "$(type -t unpad)" != function ] ; then
 unpad() {
  local string="$@"
  while [ "${string:0:1}" = " " ] ; do
   string="${string:1}"
  done
  while [ "${string:$((${#string}-1)):1}" = " " ] ; do
   string="${string:0:$((${#string}-1))}"
  done
  echo "$string"
 }
fi
if [ "$(type -t errwarn)" != function ] ; then
 errwarn() { echo "WARNING: $1" ; }
fi
if [ "$(type -t uniq_list)" != function ] ; then
 uniq_list() {
  local x ; unpad $({ for x in "$@"; do echo $x ; done ; } | sort | uniq)
 }
fi
if [ "$(type -t uncap)" != function ] ; then
 case "${BASH_VERSION%%.*}" in
 0|1|2|3)
  uncap() {
   # Turn upper case into lower case in $1.
   # NB, this is a lot faster than calling 'tr' due to the call overhead.
   local i string string_out="" n c
   i=0 ; string="$1" ; n=${#string}
   while ((i<n)) ; do c="${string:$i:1}" ; i=$((i+1))
    case "$c" in
    A) c=a ;; B) c=b ;; C) c=c ;; D) c=d ;; E) c=e ;; F) c=f ;; G) c=g ;;
    H) c=h ;; I) c=i ;; J) c=j ;; K) c=k ;; L) c=l ;; M) c=m ;; N) c=n ;;
    O) c=o ;; P) c=p ;; Q) c=q ;; R) c=r ;; S) c=s ;; T) c=t ;; U) c=u ;;
    V) c=v ;; W) c=w ;; X) c=x ;; Y) c=y ;; Z) c=z ;;
    esac
    string_out="$string_out$c"
   done
   echo "$string_out"
  } ;;
 *) uncap() { echo "${1,,}" ; } ;;
 esac
fi

# START FUNCTIONS

define_tags_init() {
 taglist_scalar="RELPATHNAMES"
 taglist_block=""
 taglist_make=""
 variable_list=""
 meta_variable_list=""
}

define_tags() {
 # Define tags:
 # Scalar tags
 taglist_scalar="\
  DESCRIPTION MAINTAINER DATE COMMENT QUEUEING_SYSTEM\
  HOSTNAME DOMAIN ARCH KERNEL DISTRIBUTION OS\
  F90_VERSION CC_VERSION CXX_VERSION\
  RUN_SINGLE RUN_PARALLEL CLUSTER_RUN_MODE RUN_CLUSTER SUBMIT_SCRIPT\
  FORCE_PATH\
  CORES_PER_NODE CORES_PER_NODE_CLUSTER\
  MIN_NCORE MAX_NCORE ALLOWED_NCORE\
  MIN_NNODE MAX_NNODE ALLOWED_NNODE\
  MIN_NNODE_ENSEMBLE\
  TIME_FORMAT\
  RELPATHNAMES\
  SCRIPTCSH\
  MAKE_EXECUTABLE\
  MIN_WALLTIME MAX_WALLTIME ALLOWED_WALLTIME WALLTIME_CODES\
  MIN_CORETIME MAX_CORETIME\
  MAX_NJOBS"
 # Block tags
 taglist_block="\
  SCRIPT_HEAD SCRIPT_RUN\
  COMMAND_CHECK_F90_VERSION COMMAND_CHECK_F90\
  COMMAND_CHECK_CC_VERSION COMMAND_CHECK_CC\
  COMMAND_CHECK_CXX_VERSION COMMAND_CHECK_CXX "
 # Makefile tags
 taglist_make="\
  TYPE SHELL MPI_VERSION\
  NEED_ETIME CFLAGS_ETIME MODNAME_BUG IS_CYGWIN NATIVE_WINDOWS\
  LIB_PATH INCLUDE_DIR\
  F90 FFLAGS_opt FFLAGS_debug FFLAGS_dev FFLAGS_prof FFLAGS_all FFLAGS_NOOPT\
  LDF90 LDFLAGS_opt LDFLAGS_debug LDFLAGS_dev LDFLAGS_prof LDFLAGS_all\
  LDLIBS_opt LDLIBS_debug LDLIBS_dev LDLIBS_prof LDLIBS_all\
  CC CFLAGS_opt CFLAGS_debug CFLAGS_dev CFLAGS_prof CFLAGS_all\
  LDCC LDCFLAGS_opt LDCFLAGS_debug LDCFLAGS_dev LDCFLAGS_prof LDCFLAGS_all\
  LDCLIBS_opt LDCLIBS_debug LDCLIBS_dev LDCLIBS_prof LDCLIBS_all\
  CXX CXXFLAGS_opt CXXFLAGS_debug CXXFLAGS_dev CXXFLAGS_prof CXXFLAGS_all\
  LDCXX LDCXXFLAGS_opt LDCXXFLAGS_debug LDCXXFLAGS_dev LDCXXFLAGS_prof\
  LDCXXFLAGS_all\
  LDCXXLIBS_opt LDCXXLIBS_debug LDCXXLIBS_dev LDCXXLIBS_prof LDCXXLIBS_all\
  CFLAGS_F90_INTERFACE AR\
  SUPPORT_OPENMP FFLAGS_OPENMP_yes FFLAGS_OPENMP_no LDFLAGS_OPENMP_yes\
  LDFLAGS_OPENMP_no CFLAGS_OPENMP_yes CFLAGS_OPENMP_no\
  SUPPORT_SHM CFLAGS_SHM\
  HAVE_BLAS HAVE_LAPACK LDBLAS_yes LDLAPACK_yes FFLAGS_libs FFLAGS0_libs\
  INSTALL_DIR ENVIRONMENT_COMMAND\
  UTILS_MODE\
  F90_NATIVE FFLAGS_opt_NATIVE FFLAGS_all_NATIVE\
   FFLAGS_libs_NATIVE FFLAGS0_libs_NATIVE\
  LDF90_NATIVE LDFLAGS_opt_NATIVE LDFLAGS_all_NATIVE\
  CC_NATIVE CFLAGS_opt_NATIVE CFLAGS_all_NATIVE CFLAGS_F90_INTERFACE_NATIVE\
  LDC_NATIVE LDCFLAGS_opt_NATIVE LDCFLAGS_all_NATIVE\
  HAVE_BLAS_NATIVE LDBLAS_yes_NATIVE LDBLAS_no_NATIVE\
  HAVE_LAPACK_NATIVE LDLAPACK_yes_NATIVE LDLAPACK_no_NATIVE\
  LDLIBS_opt_NATIVE LDLIBS_all_NATIVE\
  ENVIRONMENT_COMMAND_NATIVE\
  MATH_LIB_PATH"
 # Define list of substitutable variables
 variable_list="TYPE OUT NPROC NNODE NCORE PPN TPP TPN WALLTIME BINARY SCRIPT\
  F90 CC CXX NJOB NPROC_TOTAL NNODE_TOTAL NCORE_TOTAL BINARY_ARGS"
 meta_variable_list=""
}

clear_tags() {
 # Unset the values of all tags in list $1.
 local taglist="$1" tag subtag
 if [ ! -z "$taglist" ] ; then
  for tag in $taglist ; do
   unset $tag ${tag}_deps ${tag}_is_command
  done
  unset ERROR
 else
  taglist="$taglist_scalar $taglist_block $taglist_make"
  for tag in $taglist ; do
   unset $tag ${tag}_deps ${tag}_is_command
  done
  for tag in $user_taglist ; do # (custom_taglist from previous load)
   unset user_$tag
   for subtag in description default allowed min max ; do
    unset user_${subtag}_$tag user_${subtag}_${tag}_deps\
     user_${subtag}_${tag}_is_command
   done
  done
  unset user_taglist
  for tag in $internal_taglist ; do # (custom_taglist from previous load)
   unset internal_$tag internal_${tag}_deps internal_${tag}_is_command
  done
  unset internal_taglist
  unset ERROR
 fi
}

load_tags() {
 # Load tags from file $1.
 local file="$1" line tagline IFS_save tag ivector value is_in_list iline
 local is_error tagname deps cmd_type t1 t2 line2
 if [ ! -f "$file" ] ; then
  ERROR="$ERROR file_not_found:$file"
  return
 fi
 # Process includes first.
 { while read line ; do load_tags "${line#*/}" ; done ; }\
  < <(grep -E "^ *include " "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 # Process current file: tags.
 tag=""
 ivector=-1
 is_error=0
 is_in_list=0
 IFS_save="$IFS"
 {
  while : ; do
   IFS=""
   read -r line || break
   IFS="$IFS_save"
   iline=${line%%:*}
   line="${line#*:}"
   tagline="${line#\#-!}"
   if [ "${tagline:0:2}" = "  " ] ; then
    # Indentation for block content or continuation lines
    value=$(unpad "$tagline")
    if [ -z "$tag" ] ; then
     ((is_error==0)) && ERROR="$ERROR syntax_error:$file:$iline"
     is_error=1
     continue
    fi
    if ((ivector<0)) ; then
     deps=$(variables_in_line "$value")
     eval "$tag=\"\$$tag \$value\" ; ${tag}_deps=\"\$${tag}_deps $deps\""
    else
     ivector=$((ivector+1))
     deps=$(variables_in_line "$value")
     eval "$tag[$ivector]=\"\$value\" ; $tag[0]=$ivector ;
      ${tag}_deps=\"\$${tag}_deps $deps\""
    fi
   else
    # Declaration of new tag (indentation checked below)
    ivector=-1 ; tag="" ; is_error=0
    if [ "${tagline:0:1}" != " " ] ; then
     ERROR="$ERROR missing_space:$file:$iline"
     is_error=1
     continue
    fi
    if [ "${tagline//:/}" = "$tagline" ] ; then
     ERROR="$ERROR missing_colon:$file:$iline"
     is_error=1
     continue
    fi
    tagname="$(unpad ${tagline%%:*})"
    if [ "${tagname:0:1}" = "*" ] ; then
     tagname="${tagname:1}"
     case "$tagname" in
     USER.*)
      tagname="${tagname#*.}" ; t1="" ; t2=""
      case "$tagname" in
      DESCRIPTION.*) t1=user_description ; t2="${tagname#*.}" ;;
      DEFAULT.*) t1=user_default ; t2="${tagname#*.}" ;;
      ALLOWED.*) t1=user_allowed ; t2="${tagname#*.}" ;;
      MIN.*) t1=user_min ; t2="${tagname#*.}" ;;
      MAX.*) t1=user_max ; t2="${tagname#*.}" ;;
      esac
      if [ ! -z "$t1" ] ; then
       user_varlist="$user_varlist $t2"
       tagname="${t1}_$t2"
       cmd_type=1
      else
       cmd_type=0
      fi ;;
     INTERNAL.*)
      t1=internal ; t2="${tagname#*.}"
      internal_varlist="$internal_varlist $t2"
      tagname="${t1}_$t2"
      cmd_type=1 ;;
     *)
      if in_line $tagname $taglist_scalar ; then
       cmd_type=1
      elif in_line $tagname $taglist_block ; then
       cmd_type=2
      else
       cmd_type=0
      fi ;;
     esac
     if ((cmd_type>0)) ; then
      value="$(unpad ${tagline#*:})"
      if [ -z "$value" ] ; then
       ivector=0 ; tag=$tagname ; unset $tag ${tag}_deps
       eval "${tag}_is_command=$cmd_type"
      else
       ERROR="$ERROR bad_cmd:$file:$iline:$tagname"
       is_error=1
       continue
      fi
     else
      ERROR="$ERROR unknown_tag:$file:$iline:$tagname"
      is_error=1
      continue
     fi
    else
     case "$tagname" in
     USER.*)
      tagname="${tagname#*.}" ; t1="" ; t2=""
      case "$tagname" in
      DESCRIPTION.*) t1=user_description ; t2="${tagname#*.}" ;;
      DEFAULT.*) t1=user_default ; t2="${tagname#*.}" ;;
      ALLOWED.*) t1=user_allowed ; t2="${tagname#*.}" ;;
      MIN.*) t1=user_min ; t2="${tagname#*.}" ;;
      MAX.*) t1=user_max ; t2="${tagname#*.}" ;;
      esac
      if [ ! -z "$t1" ] ; then
       user_varlist="$user_varlist $t2"
       tagname="${t1}_$t2"
       tag=$tagname
       value="$(unpad ${tagline#*:})"
       deps=$(variables_in_line "$value")
       unset ${tag}_is_command
       eval "$tag=\"\$value\" ; ${tag}_deps=\"$deps\""
      else
       ERROR="$ERROR unknown_tag:$file:$iline:$tagname"
       is_error=1
       continue
      fi ;;
     INTERNAL.*)
      t1=internal ; t2="${tagname#*.}"
      internal_varlist="$internal_varlist $t2"
      tagname="${t1}_$t2"
      tag=$tagname
      value="$(unpad ${tagline#*:})"
      deps=$(variables_in_line "$value")
      unset ${tag}_is_command
      eval "$tag=\"\$value\" ; ${tag}_deps=\"$deps\"" ;;
     *)
      if in_line $tagname $taglist_scalar ; then
       tag=$tagname
       value="$(unpad ${tagline#*:})"
       deps=$(variables_in_line "$value")
       unset ${tag}_is_command
       eval "$tag=\"\$value\" ; ${tag}_deps=\"$deps\""
      elif in_line $tagname $taglist_block ; then
       tag=$tagname
       value="$(unpad ${tagline#*:})"
       if [ -z "$value" ] ; then
        ivector=0
        unset $tag ${tag}_is_command ${tag}_deps
       else
        ERROR="$ERROR bad_block:$file:$iline:$tagname"
        is_error=1
        continue
       fi
      else
       ERROR="$ERROR unknown_tag:$file:$iline:$tagname"
       is_error=1
       continue
      fi ;;
     esac
    fi
   fi
  done
 } < <(grep -nE "^#-!" "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 IFS="$IFS_save"
 # Process current file: make tags.
 {
  while : ; do
   IFS=""
   read -r line || break
   IFS="$IFS_save"
   iline=${line%%:*}
   line="${line#*:}"
   while [ "${line:$((${#line}-1)):1}" = \\ ] ; do
    IFS=""
    read -r line2 || break
    IFS="$IFS_save"
    line2="${line2#*:}"
    line="${line:0:$((${#line}-1))} $line2"
   done
   IFS="$IFS_save"
   tag="${line%%=*}"
   [ "$tag" = "$line" ] && continue # no "=" sign?
   tag=$(unpad "$tag")
   [ -z "$tag" ] && continue # no tag name?
   if in_line $tag $taglist_make ; then
    value=$(unpad "${line#*=}")
    eval "$tag=\"\$value\""
   else
    ERROR="$ERROR unknown_maketag:$file:$iline:$tag"
   fi
  done
 } < <(grep -nvE "^ *#|^ *include |^ *$" "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 IFS="$IFS_save"
}

quickload_tags_scalar() {
 # Load scalar tags in list $1 from file $2.
 local taglist="$1" file="$2" line line2 IFS_save value deps tag cont lcont
 if [ ! -f "$file" ] ; then
  ERROR="$ERROR file_not_found:$file"
  return
 fi
 # Process includes first.
 { while read line ; do quickload_tags_scalar "$taglist" "${line#*/}" ; done\
  ; } < <(grep -E "^ *include " "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 # Process current file.
 cont="#-!  " ; lcont=${#cont}
 line2=""
 IFS_save="$IFS"
 {
  while : ; do
   if [ -z "$line2" ] ; then
    IFS=""
    read -r line || break
    IFS="$IFS_save"
    line="${line:4}"
   else
    line="${line2:4}" ; line2=""
   fi
   [ "${line:0:1}" = " " ] && continue
   tag="${line%%:*}"
   in_line "$tag" $taglist || continue
   value="$(unpad "${line#*:}")"
   while : ; do
    IFS=""
    read -r line2 || break
    IFS="$IFS_save"
    [ "${line2:0:$lcont}" = "$cont" ] || break
    value="$value $(unpad "${line2:$lcont}")"
    line2=""
   done
   IFS="$IFS_save"
   deps=$(variables_in_line "$value")
   unset ${tag}_is_command
   eval "$tag=\"\$value\" ; ${tag}_deps=\"$deps\""
  done
 } < <(grep -E "^#-!" "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 IFS="$IFS_save"
}

quickload_tags_make() {
 # Load makefile tags in list $1 from file $2.
 local taglist="$1" file="$2" line line2 IFS_save value tag
 if [ ! -f "$file" ] ; then
  ERROR="$ERROR file_not_found:$file"
  return
 fi
 # Process includes first.
 { while read line ; do quickload_tags_make "$taglist" "${line#*/}" ; done ; }\
  < <(grep -E "^ *include " "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 lhead=${#tag}
 IFS_save="$IFS"
 {
  while : ; do
   IFS=""
   read -r line || break
   IFS="$IFS_save"
   while [ "${line:$((${#line}-1)):1}" = \\ ] ; do
    IFS=""
    read -r line2 || break
    IFS="$IFS_save"
    line="${line:0:$((${#line}-1))} $line2"
   done
   IFS="$IFS_save"
   tag="$(unpad "${line%%=*}")"
   in_line "$tag" $taglist || continue
   value=$(unpad "${line#*=}")
   eval "$tag=\"\$value\""
  done
 } < <(grep -vE "^ *#|^ *include |^ *$" "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 IFS="$IFS_save"
}

quickload_tags_block() {
 # Load block tags in list $1 from file $2.
 local taglist="$1" file="$2" line line2 IFS_save value deps tag cont lcont
 local ivector
 if [ ! -f "$file" ] ; then
  ERROR="$ERROR file_not_found:$file"
  return
 fi
 # Process includes first.
 { while read line ; do quickload_tags_block "$taglist" "${line#*/}" ; done\
  ; } < <(grep -E "^ *include " "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 # Process current file.
 cont="#-!  " ; lcont=${#cont}
 IFS_save="$IFS"
 {
  while : ; do
   if [ -z "$line2" ] ; then
    IFS=""
    read -r line || break
    IFS="$IFS_save"
    line="${line:4}"
   else
    line="${line2:4}" ; line2=""
   fi
   [ "${line:0:1}" = " " ] && continue
   tag="${line%%:*}"
   in_line "$tag" $taglist || continue
   unset $tag ${tag}_deps ${tag}_is_command
   ivector=0 ; eval "$tag[0]=$ivector"
   while : ; do
    IFS=""
    read -r line2 || break
    IFS="$IFS_save"
    [ "${line2:0:$lcont}" = "$cont" ] || break
    value="$(unpad "${line2:$lcont}")"
    line2=""
    ivector=$((ivector+1))
    deps=$(variables_in_line "$value")
    eval "$tag[$ivector]=\"\$value\" ; $tag[0]=$ivector ;
     ${tag}_deps=\"\$${tag}_deps $deps\""
   done
  done
 } < <(grep -E "^#-!" "$file")
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
 IFS="$IFS_save"
}

substitute_variables() {
 # Substitute variables in "$1" with their values.
 local line="$1" line_out="" line2 var
 while : ; do
  line2="${line#*&}" ; var="${line2%%&*}"
  if [ "$line" = "$line2" ] || [ "$line2" = "$var" ] ; then
   line_out="$line_out$line"
   break
  fi
  if [[ "$var" == +([A-Z0-9_.]) ]] ; then
   case "$var" in
   USER.*) eval "value=\"\$user_${var#*.}\"" ;;
   INTERNAL.*) eval "value=\"\$internal_${var#*.}\"" ;;
   ENV.*) eval value="\$${var#*.}" ;;
   *) eval "value=\"\$var_$var\"" ;;
   esac
   line_out="$line_out${line%%&$var&*}$value"
   line="${line2#*&}"
  else
   line_out="$line_out${line%%&*}&"
   line="$line2"
  fi
 done
 echo "$line_out"
}

variables_in_line() {
 # Output the names (without ampersands) of the variables found in line $1.
 local line="$1" line2 var varlist=""
 while : ; do
  line2="${line#*&}" ; var="${line2%%&*}"
  [ "$line" = "$line2" ] || [ "$var" = "$line2" ] && break
  if [[ "$var" == +([A-Z0-9_.]) ]] ; then
   varlist="$varlist $var"
   line="${line2#*&}"
  else
   line="$line2"
  fi
 done
 unpad "$varlist"
}

eval_tag() {
 # Evaluate value of tag "$1" using variable subsitution/command substitution.
 local tag="$1" value v1 nlines IFS_save="$IFS" scalar_target
 eval "value=\"\$$tag\" ; is_command=\"\$${tag}_is_command\""
 [ -z "$value" ] && return
 case "$is_command" in
 1|2) # A command
  # Evaluate variables first
  nlines=$value
  iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
   eval "$tag[$iline]=\"\$(substitute_variables \"\${$tag[$iline]}\")\""
  done
  # Now extract final tag value
  case "$is_command" in
  1) # A command that evaluates to a scalar
   {
    unset $tag ${tag}_is_command ${tag}_deps
    IFS=""
    if ! read -r line ; then
     IFS="$IFS_save"
     return
    fi
    IFS="$IFS_save"
    eval "$tag=\"\$(unpad \"\$line\")\""
    while : ; do
     IFS=""
     read -r line || break
     IFS="$IFS_save"
     eval "$tag=\"\$$tag \$(unpad \"\$line\")\""
    done
   } < <(eval_code $tag)
   (exit 0) # bypass bash fd leak (v3.2 - v4.1)
   IFS="$IFS_save" ;;
  2) # A command that evaluates to a block
   nlines=0
   {
    unset $tag ${tag}_is_command ${tag}_deps
    while : ; do
     IFS=""
     read -r line || break
     IFS="$IFS_save"
     nlines=$((nlines+1))
     eval "$tag[$nlines]=\"\$line\" ; $tag[0]=$nlines"
    done
   } < <(eval_code $tag)
   (exit 0) # bypass bash fd leak (v3.2 - v4.1)
   IFS="$IFS_save" ;;
  esac ;;
 *) # Not a command
  eval "v1=\"\${$tag[1]}\""
  if [ -z "$v1" ] ; then # scalar
   eval "$tag=\"\$(substitute_variables \"\$value\")\""
  else # block
   nlines=$value
   iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
    eval "$tag[$iline]=\"\$(substitute_variables \"\${$tag[$iline]}\")\""
   done
  fi ;;
 esac
}

print_block_tag() {
 # Print the contents in block tag $1 to stdout.
 # Use eval_tag first to perform variable/command substitution.
 local tag="$1" nlines iline line
 eval "nlines=\"\${$tag[0]}\""
 [ -z "$nlines" ] && return
 iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
  eval "echo \"\${$tag[$iline]}\""
 done
}

print_block_tag_csh() {
 # Print the contents in block tag $1 to stdout.
 # Use eval_tag first to perform variable/command substitution.
 # Adds a background ampersand required for csh batch scripts.
 local tag="$1" nlines iline line
 eval "nlines=\"\${$tag[0]}\""
 [ -z "$nlines" ] && return
 iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
  eval "echo \"\${$tag[$iline]}\" \&"
 done
}

eval_code() {
 # Evaluate the code in block tag $1, passing the output on standard output.
 # Use eval_tag first to perform variable/command substitution.
 local tag="$1" nlines iline line
 eval nlines=\${$tag[0]}
 [ -z "$nlines" ] && return
 bash < <(
  iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
   eval "echo \"\${$tag[$iline]}\""
  done
 ) 2> /dev/null
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
}

eval_code_tmp() {
 # Create and change into a temporary directory, evaluate the code in block
 # tag $1, and delete the temporary directory, passing the output on standard
 # output. Use eval_tag first to perform variable/command substitution.
 local tag="$1" nlines iline line
 eval nlines=\${$tag[0]}
 [ -z "$nlines" ] && return
 bash < <(
  echo "eval_code_tmp_dir=\"\$(mktemp -d 2>/dev/null)\""
  echo "cd \"\$eval_code_tmp_dir\" >& /dev/null"
  iline=0 ; while ((iline<nlines)) ; do iline=$((iline+1))
   eval "echo \"\${$tag[$iline]}\""
  done
  echo "cd >& /dev/null"
  echo "rm -rf \"\$eval_code_tmp_dir\" >& /dev/null"
 ) 2> /dev/null
 (exit 0) # bypass bash fd leak (v3.2 - v4.1)
}

report_errors() {
 # Report errors.
 local errcode errfile errline errkey missing_tags tag
 local req_taglist fvalue nlines
 if [ ! -z "$ERROR" ] ; then
  for error in $ERROR ; do
   errcode="" ; errfile="" ; errline="" ; errkey=""
   errcode=${error%%:*} ; errtail=${error#*:}
   if [ "$error" != "$errtail" ] ; then
    error="$errtail" ; errfile=${error%%:*} ; errtail=${error#*:}
    if [ "$error" != "$errtail" ] ; then
     error="$errtail" ; errline=${error%%:*} ; error=${error#*:}
     [ "$error" != "$errtail" ] && errkey=$error
    fi
   fi
   case "$errcode" in
   file_not_found) errwarn "file $errfile not found" ;;
   syntax_error) errwarn "syntax error at line $errline of $errfile" ;;
   missing_space) errwarn "missing space after '#-!' at line $errline of\
    $errfile" ;;
   missing_colon) errwarn "missing colon at line $errline of $errfile" ;;
   bad_block) errwarn "declaration of block tag '$errkey' at line $errline of\
    $errfile contains non-blanks after the colon" ;;
   bad_cmd) errwarn "declaration of command for tag '$errkey' at line $errline\
    of $errfile contains non-blanks after the colon" ;;
   unknown_tag) errwarn "unknown tag '$errkey' at line $errline of\
    $errfile" ;;
   unknown_maketag) errwarn "unknown Makefile variable '$errkey' at line\
    $errline of $errfile" ;;
   *) errwarn "unhandled error found identified by string '$error'" ;;
   esac
  done
 fi
 for req in mandatory recommended ; do
  eval "req_taglist=\"\$${req}_taglist\""
  missing_tags=""
  for tag in $req_taglist ; do
   eval "fvalue=\"\$$tag\""
   [ -z "$fvalue" ] && missing_tags="$missing_tags $tag"
  done
  if [ ! -z "$missing_tags" ] ; then
   if (($(set -- $missing_tags ; echo $#)==1)) ; then
    errwarn "$req tag '$(unpad $missing_tags)' missing for\
     CASINO_ARCH=$CASINO_ARCH."
   else
    errwarn "$req tags '$(unpad $missing_tags)' missing for\
     CASINO_ARCH=$CASINO_ARCH."
   fi
  fi
 done
}

merge_vardeps() {
 # Merge a set of variables $2 in the dependency tree so that they are all
 # evaluated together under the name META.$1.
 local meta_name="$1" varlist="$2" var
 meta_variable_list="$meta_variable_list $meta_name"
 eval "meta_vars_$meta_name=\"$varlist\""
 for var in $varlist ; do
  eval "map_to_meta_$var=$meta_name ; need_$var=0"
 done
}

declare_vardeps() {
 # Declare the tags $2 on which variable/metavariable $1 depends, and the
 # command to evaluate it $3 (usually a function name).
 local var="$1" taglist="$2" function_eval="$3"
 case "$var" in
 META.*) var=meta_${var#*.} ;;
 esac
 eval "vardepends_$var=\"\$taglist\" ;\
  vareval_$var=\"\$function_eval\""
}

resolve_dependencies() {
 # Determine which tags depend on which variables, which variables depend
 # on which tags, establish in which order they must be evaluated in
 # order to end up evaluating the tags in list $1, and evaluate everything.
 local target_tags="$1" default allowed min max ilevel jlevel tag var taglist
 local varlist deps varlist2 taglist2 t1 t2 nlevel map infocmd="$2"
 local -a tags_in_level vars_in_level
 [ -z "$infocmd" ] && infocmd=:
 tags_in_level[1]="$target_tags"
 ilevel=1 ; while : ; do
  # Get list of dependencies of current tags
  varlist=""
  for tag in ${tags_in_level[$ilevel]} ; do
   case "$tag" in
   user_*|internal_*) : ;;
   *) in_line $tag $taglist_scalar $taglist_block || { errwarn "Encountered\
    unknown tag '$tag' when resolving dependencies." ; return 1 ; } ;;
   esac
   eval "varlist2=\"\$${tag}_deps\""
   for var in $varlist2 ; do
    case "$var" in
    *.*) : ;;
    *)
     eval "map=\"\$map_to_meta_$var\""
     if [ ! -z "$map" ] ; then
      eval "need_$var=1"
      var=META.$map
     fi ;;
    esac
    varlist="$varlist $var"
   done
  done
  vars_in_level[$ilevel]="$varlist"
  [ -z "$(unpad $varlist)" ] && break
  # "Float" variables to current level if they are repeated
  jlevel=0 ; while ((jlevel<ilevel-1)) ; do jlevel=$((jlevel+1))
   varlist2="${vars_in_level[$jlevel]}"
   for var in $varlist ; do
    in_line $var $varlist2 && varlist2="$(rem_list $var $varlist2)"
   done
   [ -z "$varlist2" ] && { errwarn "Cyclic dependencies between\
    tags/variables." ; return 1 ; }
   vars_in_level[$jlevel]="$varlist2"
  done
  ilevel=$((ilevel+1))
  # Get list of dependencies of current variables
  taglist=""
  for var in ${vars_in_level[$((ilevel-1))]} ; do
   case "$var" in
   USER.*)
    eval value="\$user_${var#*.}"
    taglist2="user_allowed_${var#*.} user_min_${var#*.}\
     user_max_${var#*.}"
    [ -z "$value" ] && taglist2="$taglist2 user_default_${var#*.}" ;;
   INTERNAL.*) taglist2="internal_${var#*.}" ;;
   ENV.*) taglist2="" ;;
   META.*) eval "taglist2=\"\$vardepends_meta_${var#*.}\"" ;;
   *)
    in_line $var $variable_list || { errwarn "Encountered unknown variable\
     '$var' when resolving dependencies." ; return 1 ; }
    eval "taglist2=\"\$vardepends_$var\"" ;;
   esac
   taglist="$taglist $taglist2"
  done
  tags_in_level[$ilevel]="$taglist"
  [ -z "$(unpad $taglist)" ] && break
  # "Float" tags to current level if they are repeated (except from level 1)
  jlevel=1 ; while ((jlevel<ilevel-1)) ; do jlevel=$((jlevel+1))
   taglist2="${tags_in_level[$jlevel]}"
   for tag in $taglist ; do
    in_line $tag $taglist2 && taglist2="$(rem_list $tag $taglist2)"
   done
   [ -z "$taglist2" ] && { errwarn "Cyclic dependencies between\
    tags/variables." ; return 1 ; }
   tags_in_level[$jlevel]="$taglist2"
  done
 done
 nlevel=$ilevel
 # Make lists unique
 $infocmd "Dependency tree:"
 ilevel=$((nlevel+1)) ; while ((ilevel>1)) ; do ilevel=$((ilevel-1))
  varlist="$(uniq_list ${vars_in_level[$ilevel]})"
  vars_in_level[$ilevel]="$varlist"
  [ -z "$varlist" ] || $infocmd "vars[$ilevel] = $varlist"
  taglist="$(uniq_list ${tags_in_level[$ilevel]})"
  tags_in_level[$ilevel]="$taglist"
  [ -z "$taglist" ] || $infocmd "tags[$ilevel] = $taglist"
 done
 # Evaluate everything
 ilevel=$((nlevel+1)) ; while ((ilevel>1)) ; do ilevel=$((ilevel-1))
  # Evaluate variables
  for var in ${vars_in_level[$ilevel]} ; do
   value=""
   case "$var" in
   USER.*)
    var="${var#*.}"
    eval "default=\"\$user_default_$var\""
    eval "allowed=\"\$user_allowed_$var\""
    eval "min=\"\$user_min_$var\""
    check_number_N "$min" || min=""
    eval "max=\"\$user_max_$var\""
    check_number_N "$max" || max=""
    eval "value=\"\$user_$var\""
    [ -z "$value" ] && [ ! -z "$default" ] && value="$default"
    if [ -z "$value" ] ; then
     case "$min.$max.$allowed" in
     ..) : ;;
     ..*) value="$(field 1 $allowed)" ;;
     *.*.) [ ! -z "$min" ] && value="$min" || value="$max" ;;
     *.*.*)
      for t1 in $allowed ; do
       check_number_N "$t1" || continue
       [ ! -z "$min" ] && ((t1<min)) && continue
       [ ! -z "$max" ] && ((t1>max)) && continue
       [ -z "$value" ] || ((t1<value)) && value="$t1"
      done ;;
     esac
    fi
    # Check allowed values and min/max
    [ -z "$value" ] && { errwarn "Could not assign a value to USER.$var\
     from the default/allowed/min/max values in $CASINO_ARCH.arch. Provide\
     its value on the command line with --user.$(uncap $var)=<value>.\
     Information on possible values may be available by typing 'runqmc --help'.\
     " ;\
     return 1 ; }
    [ ! -z "$allowed" ] && ! in_line "$value" $allowed && { errwarn "Value of\
     variable USER.$var='$value' not in allowed-value list '$allowed'." ;\
     return 1 ; }
    [ ! -z "$min$max" ] && ! check_number_N "$value" && { errwarn "Value of\
    variable USER.$var='$value' is not a positive integer, which it is\
    required to be since a minimum and/or a maximum value has been defined\
    in $CASINO_ARCH.arch" ; return 1 ; }
    [ ! -z "$min" ] && ((value<min)) && { errwarn "Value of variable\
     USER.$var='$value' is below minimum of $min." ; return 1 ; }
    [ ! -z "$max" ] && ((value>max)) && { errwarn "Value of variable\
     USER.$var='$value' is above maximum of $max." ; return 1 ; }
    eval "user_$var=\"$value\"" ;;
   INTERNAL.*) eval "value=\"\$internal_${var#*.}\"" ;;
   ENV.*) eval "value=\"\$${var#*.}\"" ;;
   META.*)
    var="meta_${var#*.}"
    eval "vareval=\"\$vareval_$var\""
    [ -z "$vareval" ] || eval "$vareval" ;;
   *)
    eval "vareval=\"\$vareval_$var\""
    [ -z "$vareval" ] || eval "$vareval" ;;
   esac
   [ -z "$value" ] || $infocmd "Evaluated $var='$value'"
  done
  # Evaluate tags
  for tag in ${tags_in_level[$ilevel]} ; do
   eval_tag $tag
   if in_line $tag $taglist_block ; then
    $infocmd "Evaluated $tag:"
    ((verbosity>=2)) && print_block_tag $tag
   else
    eval value="\$$tag"
    $infocmd "Evaluated $tag='$value'"
   fi
  done
 done
 return 0
}
