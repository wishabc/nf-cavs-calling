#!/bin/bash
# Usage: extract_symlink <target-dir> <-f | -n>
# -f moves symlink targets safely
# -n performs a 'dry run', printing all the commands to be performed

set -e  # Exit on error

function extract_symlink () {
    local link="$1"
    
    if [ -h "$link" ]; then
        local target
        target=$(realpath "$link")

        if [ -e "$target" ]; then
            local tmpdir
            tmpdir=$(mktemp -d --suffix="_symlink_extract")

            cp -r "$target" "$tmpdir/" && rm "$link" && mv "$tmpdir/"* "$link"
            rmdir "$tmpdir"  # Remove temp directory after use
        else
            echo "Error: Target '$target' does not exist" >&2
        fi
    fi
}

export -f extract_symlink

if [ ! -d "$1" ]; then
    echo "Error: Target directory '$1' does not exist or is not a directory" >&2
    exit 1
fi

case "$2" in
    "-n")
        find "$1" -type l -print0 | while IFS= read -r -d '' link; do
            target=$(realpath "$link")
            if [[ -e "$target" ]]; then
                dir=$(dirname "$link")
                echo "cp -r \"$target\" \"$dir/tmpfile\" && rm \"$link\" && mv \"$dir/tmpfile\" \"$link\""
            else
                echo "Error: Target '$target' does not exist" >&2
            fi
        done
        ;;
    "-f")
        find "$1" -type l -print0 | xargs -0 -I {} bash -c 'extract_symlink "$@"' _ {}
        ;;
    *)
        echo "Error: Neither -f nor -n are provided" >&2
        exit 1
        ;;
esac
