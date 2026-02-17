#!/bin/bash
#
# Gzips all .fits files from the .galactic group subdirectory and places the
# resulting .fits.gz files into the main group-level delivery directory.
#
# Usage: bash gzip_files.sh [--execute]
#
# Dry-run by default. Pass --execute to actually gzip files.

set -euo pipefail

BASE="/orange/adamginsburg/ACES/products_for_ALMA"
SRC_DIR="${BASE}/group.uid___A001_X1590_X30a9.lp_slongmore.galactic"
DST_DIR="${BASE}/group.uid___A001_X1590_X30a9.lp_slongmore"

DRY_RUN=true
if [[ "${1:-}" == "--execute" ]]; then
    DRY_RUN=false
else
    echo "=== DRY RUN MODE (pass --execute to gzip for real) ==="
fi

# Check directories exist
if [[ ! -d "$SRC_DIR" ]]; then
    echo "ERROR: Source directory not found: $SRC_DIR"
    exit 1
fi

if [[ ! -d "$DST_DIR" ]]; then
    echo "ERROR: Destination directory not found: $DST_DIR"
    exit 1
fi

count=0
skipped=0

for fitsfile in "$SRC_DIR"/*.fits; do
    [[ -e "$fitsfile" ]] || { echo "No .fits files found in $SRC_DIR"; exit 0; }

    basename=$(basename "$fitsfile")
    outfile="${DST_DIR}/${basename}.gz"

    # Skip if gzipped version already exists
    if [[ -f "$outfile" ]]; then
        echo "SKIP (already exists): $outfile"
        skipped=$((skipped + 1))
        continue
    fi

    if $DRY_RUN; then
        echo "WOULD GZIP: $fitsfile -> $outfile"
    else
        echo "GZIP: $fitsfile -> $outfile"
        gzip -c "$fitsfile" > "$outfile"
    fi
    count=$((count + 1))
done

echo ""
echo "Group level complete. Processed: $count | Skipped (existing): $skipped"

# --- Gzip .tt0. and .tt1. FITS files in member directories in-place ---
# Original unzipped files are first moved to a backup directory.

BACKUP_DIR="${BASE}/tt_images_unzipped"

if ! $DRY_RUN; then
    mkdir -p "$BACKUP_DIR"
fi

tt_count=0
tt_skipped=0

for memberdir in "$BASE"/member.*/; do
    [[ -d "$memberdir" ]] || continue

    for fitsfile in "$memberdir"*.fits; do
        [[ -e "$fitsfile" ]] || continue

        basename=$(basename "$fitsfile")

        # Only process files containing '.tt0.' or '.tt1.'
        [[ "$basename" == *".tt0."* || "$basename" == *".tt1."* ]] || continue

        gzfile="${fitsfile}.gz"

        # Skip if already gzipped in-place
        if [[ -f "$gzfile" ]]; then
            echo "SKIP (already gzipped): $gzfile"
            tt_skipped=$((tt_skipped + 1))
            continue
        fi

        if $DRY_RUN; then
            echo "WOULD BACKUP: $fitsfile -> $BACKUP_DIR/$basename"
            echo "WOULD GZIP:   $fitsfile -> $gzfile"
        else
            echo "BACKUP: $fitsfile -> $BACKUP_DIR/$basename"
            cp "$fitsfile" "$BACKUP_DIR/$basename"
            echo "GZIP:   $fitsfile -> $gzfile"
            gzip "$fitsfile"
        fi
        tt_count=$((tt_count + 1))
    done
done

echo ""
echo "Member level complete. Processed: $tt_count | Skipped (existing): $tt_skipped"