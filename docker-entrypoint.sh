#!/bin/sh
# Entry point of the cgdist image.
#   docker run IMAGE --schema ... --profiles ...   runs cgdist with these arguments
#   docker run IMAGE cgdist-cache pull ...         runs another program of the image
#   docker run IMAGE /bin/bash -c '...'            runs a shell: what Nextflow and
#                                                  other workflow managers do
if [ "$#" -eq 0 ]; then
    exec cgdist --help
fi
case "$1" in
    -*) exec cgdist "$@" ;;
esac
if command -v "$1" >/dev/null 2>&1; then
    exec "$@"
fi
exec cgdist "$@"
