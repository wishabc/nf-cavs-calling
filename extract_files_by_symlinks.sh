#!/bin/bash
# Usage: extract_symlink.sh <target-dir> <-f | -n> [num_jobs]
# -f replaces symlinks with targets (via move)
# -n prints the shell commands (dry run)
# num_jobs is optional (default = 1)
# Note: Progressbar writes to stderr

set -euo pipefail

# Build the command to replace a symlink with its target
function build_extract_command () {
    local link="$1"

    if [[ ! -h "$link" ]]; then
        return 1
    fi

    local target
    target=$(realpath "$link") || {
        echo "Error: Failed to resolve symlink: $link" >&2
        exit 1
    }

    if [[ ! -e "$target" ]]; then
        echo "Error: Target '$target' does not exist" >&2
        exit 1
    fi

    echo "unlink \"$link\" && mv \"$target\" \"$link\""
}

# Wrapper for actual execution
function extract_symlink () {
    local link="$1"

    local cmd
    cmd=$(build_extract_command "$link") || exit 1
    eval "$cmd"
}

export -f build_extract_command
export -f extract_symlink

# ---- Parse arguments ----
target_dir="$1"
mode_flag="$2" # -f or -n
num_jobs="${3:-1}"

if [[ ! -d "$target_dir" ]]; then
    echo "Error: '$target_dir' is not a directory" >&2
    exit 1
fi

if ! [[ "$num_jobs" =~ ^[0-9]+$ ]]; then
    echo "Error: num_jobs must be a positive integer" >&2
    exit 1
fi

# ---- Mode handling ----
case "$mode_flag" in
    -n) parallel_func="build_extract_command" ;;
    -f) parallel_func="extract_symlink" ;;
    *)
        echo "Error: Unknown mode flag '$mode_flag'. Use -f or -n." >&2
        exit 1
        ;;
esac

# ---- Run parallel ----
find "$target_dir" -type l -print0 \
    | parallel --bar -0 -j "$num_jobs" "$parallel_func" "{}"
