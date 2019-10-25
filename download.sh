#!/bin/sh
#
# Download a file and verify its checksum; repeat until
# the checksum matches or the user told us to give up.

set -e

usage() {
    cat <<EOF
Usage: $0 URL [SHA256 [FILENAME]]
EOF

    exit $1
}

if [ $# -lt 1 ]; then
    usage 1
fi
if [ $# -gt 3 ]; then
    usage 1
fi

FileURL="$1"
ExpectedChecksum="$2"
if [ $# = 3 ]; then
    Filename="$3"
else
    Filename=$(basename ${FileURL})
fi

while : ; do

    #
    # Download
    #
    if [ -e ${Filename} ] ; then
        echo "Checking for ${Filename}... found existing file, using that"
        # TODO: add a checksum check?
    elif command -v wget >/dev/null 2>&1 ; then
        echo "Checking for ${Filename}... downloading with wget"
        wget -O ${Filename} ${FileURL}
    elif command -v curl >/dev/null 2>&1 ; then
        echo "Checking for ${Filename}... downloading with curl"
        curl -L -o ${Filename} ${FileURL}
    else
        echo "Checking for ${Filename}... not found and unable to download"
        echo
        echo "Please manually download the following file:"
        echo "  ${FileURL}"
        echo "Put it into the current directory and re-run this script."
        exit 1
    fi

    #
    # Compute checksum
    #
    ActualChecksum=$(shasum -a 256 ${Filename})
    echo "Verifying SHA256 checksum of ${Filename}..."
    if command -v sha256sum >/dev/null 2>&1 ; then
        ActualChecksum=$(sha256sum ${Filename})
    elif command -v shasum >/dev/null 2>&1 ; then
        ActualChecksum=$(shasum -a 256 ${Filename})
    elif command -v openssl >/dev/null 2>&1 ; then
        ActualChecksum=$(openssl sha256 ${Filename})
    else
        ActualChecksum=skip
    fi

    #
    # Compare checksum
    #
    case ${ActualChecksum} in
       *${ExpectedChecksum}* )
           echo "   valid SHA256 checksum"
           break
           ;;
       *)  echo "   invalid SHA256 checksum, expected ${ExpectedChecksum} but got ${ActualChecksum}"
           echo "     retrying in 30 seconds..."
           sleep 30
           rm -f ${Filename}
           ;;
    esac

done;

exit 0
