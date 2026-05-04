#!/bin/bash

if [ -z "$BASH_VERSION" ]; then
    echo "Error: This script must be run in the BASH shell." >&2
    exit 1
fi

MAX_READ_LENGTH=151
MAX_MISMATCH=3
tempfile=$(mktemp -t temp-DBPZ-variScan.XXXXXXXXXXXX -u )

if [[ ${#tempfile} -le 20 ]]; then
    echo "Error: the 'mktemp' command isn't available." >&2
    exit 1
fi

SCRIPTDIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
WORKDIR=$( realpath . )
PURPLE='\033[0;35m'
NC='\033[0m' # No Color
echo -e "${PURPLE}== variScan Pipeline Starting ==${NC}"

# --- Function to convert to absolute path ---
get_abs_path() {
    # $1 is the path to convert
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

# --- 1. Check for correct number of arguments ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file1.fastq.gz> <file2.fastq.gz> <library.csv>"
    exit 1
fi

# Assign arguments to variables
R1_RAW="$1"
R2_RAW="$2"
LIB_RAW="$3"

# --- 2. Check existence and convert to absolute paths ---
# We check existence first, then convert.
for FILE in "$R1_RAW" "$R2_RAW" "$LIB_RAW"; do
    if [ ! -f "$FILE" ]; then
        echo -e "Error: File ${PURPLE}$FILE${NC} does not exist."
        exit 1
    fi
done

# Convert to absolute paths
R1="$(get_abs_path "$R1_RAW")"
R2="$(get_abs_path "$R2_RAW")"
LIB="$(get_abs_path "$LIB_RAW")"

# --- 3. Output status ---
echo -e "R1 Path:  ${PURPLE}$R1${NC}"
echo -e "R2 Path:  ${PURPLE}$R2${NC}"
echo -e "Library:  ${PURPLE}$LIB${NC}"


echo -e "\nRunning Alignment..."
for rno in 1 2
do
     mr2=
     if [[ $rno == 2 ]] ; then mr2=-modeR2; fi
     var_name="R${rno}"
     fqgz="${!var_name}"

     uniqfile=$tempfile.Tmp.uq.reads.gz
     gzip -cd $fqgz |awk 'NR%4==2' |sort -S 15G  |uniq -c  |gzip -1 -c > $uniqfile
     $SCRIPTDIR/bin/find-best-align -rlen $MAX_READ_LENGTH -R1 <( gzip -cd $uniqfile ) -lib "$LIB" $mr2 -outfile $tempfile.Tmp.R$rno.bin 
done

echo -e "Running End Matching..."
$SCRIPTDIR/bin/match-two-ends -lib "$LIB" -rlen $MAX_READ_LENGTH -maxMM $MAX_MISMATCH -binf1  $tempfile.Tmp.R1.bin  -binf2  $tempfile.Tmp.R2.bin  -R1  <( gzip -cd "$R1" ) -R2 <( gzip -cd "$R2" ) |gzip -1 -c > $tempfile.Tmp.restxt.gz

echo -e "Creating Spreadsheets..."
echo -e "${PURPLE}== Pipeline Complete ==${NC}"

rm -f $tempfile.Tmp.*
