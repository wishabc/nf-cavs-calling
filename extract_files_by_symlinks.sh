#!/bin/bash
# Usage: extract symlink <target-dir> <-f | -n>
# -f moves symlink targets
# -n perform a 'dry run', printing all the commands to be performed

# TODO: add graceful error handling

function extract_symlink () {
	if ! [ -d "$1" ]; then
        if [ -h "$1" ]; then
            mv $( realpath "$1" ) $1
        fi
	fi
}
export -f extract_symlink

case $2 in
    "-n")
        find $1 -type l -exec bash -c 'a="$@"; b=$( realpath $a ); if [[ "$b" != "$a" ]]; then echo mv $b $a; fi' bash {} \;
        ;;
    "-f")
        find $1 -type l -exec bash -c 'extract_symlink "$@"' bash {} \;
        ;;
    "*") 
        echo "Neither -f nor -n are provided"
esac
