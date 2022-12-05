#!/bin/bash
# Usage: extract symlink <target-dir> <-f | -n>
# -f moves symlink targets
# -n perform a 'dry run', printing all the commands to be performed


function extract_symlink () {
    echo $1
	if ! [ -d "$1" ]; then
        symlink_target=$( readlink "$1" )
        if [[ "$symlink_target" != "$1" ]]; then
            mv $symlink_target $1
        fi
	fi
}
export -f extract_symlink
case $2 in
    "-n")
        find $1 -exec echo "moving {}" \;
        ;;
    "-f")
        find $1 -exec bash -c 'extract_symlink "$@"' {} \;
        ;;
    "*") 
        echo "Neither -f nor -n are provided"
esac
