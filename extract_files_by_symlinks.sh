#!/bin/bash
# Usage: extract_symlink <target-dir> <-f | -n> [num_jobs]
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
            tmpdir=$(mktemp -d --tmpdir extract_XXXXXX_symlink)

            cp -rL "$target" "$tmpdir/" && rm "$link" && mv "$tmpdir/"* "$link"
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

num_jobs="${3:-1}"
if ! [[ "$num_jobs" =~ ^[0-9]+$ ]]; then
    echo "Error: Number of jobs must be a positive integer" >&2
    exit 1
fi

total_symlinks=$(find "$1" -type l | wc -l)

if [[ "$total_symlinks" -eq 0 ]]; then
    echo "No symlinks found in the directory."
    exit 0
fi

echo "Processing $total_symlinks symlinks..."


case "$2" in
    "-n")
        find "$1" -type l -print0 | while IFS= read -r -d '' link; do
            target=$(realpath "$link")
            if [[ -e "$target" ]]; then
                dir=$(dirname "$link")
                echo "cp -rL \"$target\" \"$dir/tmpfile\" && rm \"$link\" && mv \"$dir/tmpfile\" \"$link\""
            else
                echo "Error: Target '$target' does not exist" >&2
            fi
        done
        ;;
    "-f")
        find "$1" -type l -print0 | parallel --bar -0 -j "$num_jobs" extract_symlink "{}"
        ;;
    *)
        echo "Error: Neither -f nor -n are provided" >&2
        exit 1
        ;;
esac
