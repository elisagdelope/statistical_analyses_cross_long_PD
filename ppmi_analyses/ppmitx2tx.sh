#!/bin/bash

if [ ! -z $1 ]; then

    # PRINT USAGE
    if [ $1 == "-h" ]; then
       echo "Parse transcript names from PPMI data to be compatible with GENECODE transcript IDs."
       echo "USAGE: <script.sh> <1> <2> <~3>"
       echo "  1) Path to *.transcripts.sf files"
       echo "  2) Output filename"
       echo "  3) OPTIONAL: temporary folder (default:.)"
       exit 0
    fi

    if [ -z $3 ]; then
        TEMP_PATH="."
    else
        TEMP_PATH=$3
    fi

    tmp0="${TEMP_PATH}/0.tmp"
    tmp1="${TEMP_PATH}/1.tmp"
    touch ${tmp0}
    for filename in $1/*.transcripts.sf; do
        echo "Processing ${filename}"
        cut -f1 ${filename} | sed '1d' >> ${tmp0}
        sort ${tmp0} | uniq > ${tmp1}
        mv ${tmp1} ${tmp0}
    done

    echo "files tmp0 and tmp1 are $tmp0 and $tmp1"

    echo "Extracting tx and gene names"
    cut -f1 -d"|" ${tmp0} > ${tmp1}
    sed -i '' '1s/^/PPMI_TXNAME\n/' ${tmp0}
    sed -i '' '1s/^/TXNAME\n/' ${tmp1}
    paste ${tmp0} ${tmp1} > $2
    rm ${tmp0} ${tmp1}
    echo "File written to -> $2"
else
    echo "Missing input args. Use -h to print help."
fi
