#!/bin/bash
TARGET=$1
CONTROL=$2
OUTDIR=$3
GENOME=$4
NAME=$5


if [[ "$CONTROL" == "None" ]]; then
    macs2 callpeak -t $TARGET -g $GENOME -f BAMPE -n $NAME --bdg \ 
    --outdir $OUTDIR -q 0.1 --keep-dup all 2>$OUTDIR/${NAME}.summary.txt
else
    macs2 callpeak -t $TARGET -c $CONTROL -g $GENOME -f BAMPE -n $NAME --bdg \
    --outdir $OUTDIR -q 0.1 --keep-dup all 2>$OUTDIR/${NAME}.summary.txt
fi

if [[ $? -eq 0 ]]; then
    echo "MACS2 Peak Calling Success: $NAME"
else
    echo "MACS2 Peak Calling Failed: $NAME" >&2
    exit 1
fi

