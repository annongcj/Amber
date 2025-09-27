#! /bin/bash
script_dir=$(dirname $0)

source $script_dir/amber.sh

if [ -z "$SHELL" ] ; then
	bash
else
	$SHELL
fi
