# Define a couple of routines for pretty-printing in Makefiles.
# This code contains no single-quotes, so one can do
# bash -c '$(BASH_PRETTY) ; <rest-of-code>'
# This code contains no '#' symbols since some implementations of "make"
# (e.g. AIX's) fail to escape them properly. This prevents the use of
# ${#var}, so we have a "len" function to replace it (log scaling, tuned
# for strings of size 9-16).

BASH_PRETTY = \
 term_codes="" ;\
 [[ "$$TERM" == xterm* ]] && export TERM=xterm ;\
 if [ -n "$$TERM" ] && type -P tput >& /dev/null\
  && tput -S < /dev/null >& /dev/null ; then\
  term_codes="default bold black blue green cyan red purple brown grey" ;\
  {\
   IFS_save="$$IFS" ;\
   IFS=$$(echo -e "\t") ;\
   read default bold black blue green cyan red purple brown grey\
    term_columns ;\
   IFS="$$IFS_save" ;\
  } < <(echo -e "sgr0 \nht\n bold \nht\n setf 0 \nht\n\
   setf 1 \nht\n setf 2 \nht\n setf 3 \nht\n setf 4 \nht\n setf 5 \nht\n\
   setf 6 \nht\n setf 7 \n" | tput -S) ;\
 fi ;\
 len() {\
  local d=16 s="$$1" snext n=0 ;\
  [ -z "$$s" ] && echo 0 && return ;\
  while : ; do\
   snext="$${s:$$d}" ; ((n=n+d)) ;\
   [ -z "$$snext" ] && break ;\
   s="$$snext" ; ((d=d+d)) ;\
  done ;\
  ((d=d/2)) ;\
  for ((1; d>=1; d=d/2)) ; do\
   snext="$${s:$$d}" ;\
   [ -z "$$snext" ] && ((n=n-d)) || s="$$snext" ;\
  done ;\
  echo $$n ;\
 } ;\
 pretty_print_buffer_init() {\
  unset pp_line ;\
  npp_line=0 ;\
 } ;\
 pretty_print_buffer_add() {\
  npp_line=$$((npp_line+1)) ;\
  pp_line[$$npp_line]="$$1" ;\
 } ;\
 pretty_print_buffer() {\
  local ind1=$$1 ind2=$$2 lwidth=79 line="" nline=0 ;\
  local cat_codes="" nword c sind1 sind2 ;\
  sind1="$$(printf "%$${ind1}s")" ;\
  sind2="$$(printf "%$${ind2}s")" ;\
  shift 2 ;\
  while [ -n "$$1" ] ; do\
   if [ "$${1:0:2}" = "^^" ] ; then\
    eval "c=\$$$${1:2}" ;\
    cat_codes="$$cat_codes$$c" ;\
    line="$$line$$c" ;\
    shift ; continue ;\
   fi ;\
   nword=$$(len "$$1") ;\
   if ((nline==0)) ; then\
    line="$$cat_codes$$sind1$$1" ;\
    nline=$$((ind1+nword)) ;\
   else\
    if ((nline+1+nword>lwidth)) ; then\
     pretty_print_buffer_add "$$line$$default" ;\
     line="$$cat_codes$$sind2$$1" ;\
     nline=$$((ind2+nword)) ;\
    else\
     line="$$line $$1" ;\
     nline=$$((nline+1+nword)) ;\
    fi ;\
   fi ;\
   shift ;\
  done ;\
  ((nline==0)) || pretty_print_buffer_add "$$line$$default" ;\
 } ;\
 pretty_print_buffer_flush() {\
  local i ;\
  i=0 ; while ((i<npp_line)) ; do i=$$((i+1)) ;\
   echo "$${pp_line[$$i]}" ;\
  done ;\
  echo -n "$$default" ;\
  npp_line=0 ;\
 } ;\
 info_title() {\
  pretty_print_buffer_init ;\
  pretty_print_buffer 0 0 ^^brown $$1 ;\
 } ;\
 info_line() {\
  pretty_print_buffer 2 6 ^^purple "*" ^^cyan $$1 ;\
 } ;\
 info_colon() {\
  pretty_print_buffer 2 6 ^^purple "*" ^^brown $$1 ^^purple : ^^cyan $$2 ;\
 } ;\
 info_equal() {\
  pretty_print_buffer 2 6 ^^purple "*" ^^brown $$1 ^^purple = ^^cyan $$2 ;\
 } ;\
 info_end() { pretty_print_buffer_flush ; } ;\
 compile_info() {\
  local line="" ;\
  line="  $$red$$(printf "%-6s" "$$1")$$default  $$2" ;\
  shift 2;\
  while [ -n "$$1" ] ; do\
   line="$$line $$brown>$$default $$1" ;\
   shift ;\
  done ;\
  echo "$$line$$default" ;\
 } ;\
 nocompile_info() {\
  local line="" ;\
  line="  $$bold$$black$$(printf "%-6s" "skip")  $$2" ;\
  shift 2;\
  while [ -n "$$1" ] ; do\
   line="$$line > $$1" ;\
   shift ;\
  done ;\
  echo "$$line$$default" ;\
 }
