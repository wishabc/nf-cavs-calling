#!/bin/bash
# Usage: extract symlink <target-dir> <-f | -n>
# -f moves symlink targets
# -n perform a 'dry run', printing all the commands to be performed


function extract_symlink () {
	if ! [ -d "$file" ]; then
        symlink_target=$( readlink "$file" )
        if [[ "$symlink_target" != "$file" ]]; then
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
        echo Moving $1
        find $1 -exec bash -c 'extract_symlink "$0"' {} \;
        ;;
    "*") 
        echo "Neither -f nor -n are provided"
esac
