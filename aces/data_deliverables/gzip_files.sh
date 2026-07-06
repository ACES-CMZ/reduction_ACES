#!/bin/bash
#
# Gzips all .fits files using pigz and SLURM array jobs.
# Files in .galactic directory and files with 'pv' in filename are tar.gz'ed.
# Files are batched by size (max 500 GB per job).
#
# Usage: bash gzip_files.sh [--execute]
#
# Dry-run by default. Pass --execute to actually submit jobs.

set -euo pipefail

BASE="/orange/adamginsburg/ACES/products_for_ALMA"
SRC_DIR="${BASE}/group.uid___A001_X1590_X30a9.lp_slongmore.galactic"
DST_DIR="${BASE}/group.uid___A001_X1590_X30a9.lp_slongmore"
BACKUP_DIR="${BASE}/tt_images_unzipped"

# Batch size in bytes (500 GB)
MAX_BATCH_SIZE=$((500 * 1024 * 1024 * 1024))

DRY_RUN=true
if [[ "${1:-}" == "--execute" ]]; then
    DRY_RUN=false
else
    echo "=== DRY RUN MODE (pass --execute to submit jobs for real) ==="
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

# Create work directory for batch files
WORK_DIR="${BASE}/.pigz_batches_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"

# ============================================================================
# Part 1: Group-level files (compress from .galactic to parent directory)
# ============================================================================

GROUP_BATCH_LIST="${WORK_DIR}/group_batches.txt"
> "$GROUP_BATCH_LIST"

echo "=== Scanning group-level files in $SRC_DIR ==="

batch_num=0
current_batch_size=0
batch_file="${WORK_DIR}/group_batch_${batch_num}.txt"
> "$batch_file"

skipped=0
total_files=0

for fitsfile in "$SRC_DIR"/*.fits; do
    [[ -e "$fitsfile" ]] || { echo "No .fits files found in $SRC_DIR"; break; }

    basename=$(basename "$fitsfile")
    outfile="${DST_DIR}/${basename}.tgz"

    # Skip if .tgz version already exists
    if [[ -f "$outfile" ]]; then
        echo "SKIP (already exists): $outfile"
        # Remove the .fits file if .tgz exists to avoid duplicates
        if [[ -f "$fitsfile" ]]; then
            echo "  Removing duplicate .fits file: $fitsfile"
            rm "$fitsfile"
        fi
        skipped=$((skipped + 1))
        continue
    fi

    # Get file size
    filesize=$(stat -c%s "$fitsfile" 2>/dev/null || stat -f%z "$fitsfile" 2>/dev/null)
    
    # Check if we need a new batch
    if (( current_batch_size + filesize > MAX_BATCH_SIZE && current_batch_size > 0 )); then
        echo "$batch_file" >> "$GROUP_BATCH_LIST"
        batch_num=$((batch_num + 1))
        batch_file="${WORK_DIR}/group_batch_${batch_num}.txt"
        > "$batch_file"
        current_batch_size=0
        echo "Starting new group batch: $batch_num"
    fi

    # Add file to current batch
    echo "$fitsfile|$outfile" >> "$batch_file"
    current_batch_size=$((current_batch_size + filesize))
    total_files=$((total_files + 1))
done

# Add final batch if it has files
if [[ -s "$batch_file" ]]; then
    echo "$batch_file" >> "$GROUP_BATCH_LIST"
fi

GROUP_BATCHES=$(wc -l < "$GROUP_BATCH_LIST" 2>/dev/null || echo 0)
echo "Group level: $total_files files to process in $GROUP_BATCHES batches | Skipped: $skipped"

# ============================================================================
# Part 2: Member-level .tt0. and .tt1. files (compress in-place with backup)
# ============================================================================

MEMBER_BATCH_LIST="${WORK_DIR}/member_batches.txt"
> "$MEMBER_BATCH_LIST"

echo ""
echo "=== Scanning member-level .tt0. and .tt1. files ==="

batch_num=0
current_batch_size=0
batch_file="${WORK_DIR}/member_batch_${batch_num}.txt"
> "$batch_file"

tt_skipped=0
tt_total=0

for memberdir in "$BASE"/member.*/; do
    [[ -d "$memberdir" ]] || continue

    for fitsfile in "$memberdir"*.fits; do
        [[ -e "$fitsfile" ]] || continue

        basename=$(basename "$fitsfile")

        # Only process files containing '.tt0.' or '.tt1.'
        [[ "$basename" == *".tt0."* || "$basename" == *".tt1."* ]] || continue

        # Check if this is a PV file (should be .tgz'ed) - case insensitive
        if [[ "$basename" =~ [Pp][Vv] ]]; then
            gzfile="${fitsfile}.tgz"
        else
            gzfile="${fitsfile}.gz"
        fi

        # Skip if already compressed
        if [[ -f "$gzfile" ]]; then
            echo "SKIP (already compressed): $gzfile"
            # If this is a PV file and .tgz exists, remove the .fits to avoid duplicates
            if [[ "$basename" =~ [Pp][Vv] && -f "$fitsfile" ]]; then
                echo "  Removing duplicate PV .fits file: $fitsfile"
                rm "$fitsfile"
            fi
            tt_skipped=$((tt_skipped + 1))
            continue
        fi

        # Get file size
        filesize=$(stat -c%s "$fitsfile" 2>/dev/null || stat -f%z "$fitsfile" 2>/dev/null)
        
        # Check if we need a new batch
        if (( current_batch_size + filesize > MAX_BATCH_SIZE && current_batch_size > 0 )); then
            echo "$batch_file" >> "$MEMBER_BATCH_LIST"
            batch_num=$((batch_num + 1))
            batch_file="${WORK_DIR}/member_batch_${batch_num}.txt"
            > "$batch_file"
            current_batch_size=0
            echo "Starting new member batch: $batch_num"
        fi

        # Determine compression type
        if [[ "$basename" =~ [Pp][Vv] ]]; then
            comp_type="targz"
        else
            comp_type="gz"
        fi

        # Add file to current batch (format: source|backup_dest|gzfile|comp_type)
        echo "$fitsfile|$BACKUP_DIR/$basename|$gzfile|$comp_type" >> "$batch_file"
        current_batch_size=$((current_batch_size + filesize))
        tt_total=$((tt_total + 1))
    done
done

# Add final batch if it has files
if [[ -s "$batch_file" ]]; then
    echo "$batch_file" >> "$MEMBER_BATCH_LIST"
fi

MEMBER_BATCHES=$(wc -l < "$MEMBER_BATCH_LIST" 2>/dev/null || echo 0)
echo "Member level: $tt_total files to process in $MEMBER_BATCHES batches | Skipped: $tt_skipped"

# ============================================================================
# Create SLURM job script
# ============================================================================

SLURM_SCRIPT="${WORK_DIR}/pigz_worker.sh"

cat > "$SLURM_SCRIPT" << 'EOFSLURM'
#!/bin/bash
#SBATCH --job-name=pigz_fits
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/%x_%A_%a.log
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept

module load gcc/5.2.0 pigz/2.4

set -euo pipefail

# Get batch type and number from environment or parameters
BATCH_TYPE="${1}"
BATCH_FILE="${2}"
BACKUP_DIR="${3}"

echo "=== Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ==="
echo "Batch type: ${BATCH_TYPE}"
echo "Batch file: ${BATCH_FILE}"
echo "Started at: $(date)"

if [[ ! -f "$BATCH_FILE" ]]; then
    echo "ERROR: Batch file not found: $BATCH_FILE"
    exit 1
fi

processed=0
failed=0

if [[ "$BATCH_TYPE" == "group" ]]; then
    # Process group-level files (tar.gz with pigz into .tgz)
    while IFS='|' read -r fitsfile outfile; do
        echo "Processing: $fitsfile -> $outfile"
        basename=$(basename "$fitsfile")
        # Create .tgz using pigz for compression
        if tar -C "$(dirname "$fitsfile")" -cf - "$basename" | pigz -p ${SLURM_CPUS_PER_TASK} > "$outfile"; then
            # Verify the .tgz file was created and is not empty
            if [[ -s "$outfile" ]]; then
                echo "  ✓ Created: $outfile"
                # Remove the original .fits file to avoid duplicates
                if rm "$fitsfile"; then
                    echo "  ✓ Removed original: $fitsfile"
                else
                    echo "  WARNING: Failed to remove original $fitsfile"
                fi
                processed=$((processed + 1))
            else
                echo "ERROR: Created .tgz is empty or missing: $outfile"
                rm -f "$outfile"
                failed=$((failed + 1))
            fi
        else
            echo "ERROR: Failed to create .tgz for $fitsfile"
            rm -f "$outfile"
            failed=$((failed + 1))
        fi
    done < "$BATCH_FILE"
    
elif [[ "$BATCH_TYPE" == "member" ]]; then
    # Process member-level files (with backup)
    mkdir -p "$BACKUP_DIR"
    
    while IFS='|' read -r fitsfile backup_dest gzfile comp_type; do
        echo "Backing up: $fitsfile -> $backup_dest"
        if cp "$fitsfile" "$backup_dest"; then
            if [[ "$comp_type" == "targz" ]]; then
                echo "Creating .tgz in-place: $fitsfile -> $gzfile"
                basename=$(basename "$fitsfile")
                if tar -C "$(dirname "$fitsfile")" -cf - "$basename" | pigz -p ${SLURM_CPUS_PER_TASK} > "$gzfile"; then
                    # Verify the .tgz file was created and is not empty
                    if [[ -s "$gzfile" ]]; then
                        echo "  ✓ Created: $gzfile"
                        # Remove the original .fits file
                        if rm "$fitsfile"; then
                            echo "  ✓ Removed original: $fitsfile"
                        else
                            echo "  WARNING: Failed to remove original $fitsfile"
                        fi
                        processed=$((processed + 1))
                    else
                        echo "ERROR: Created .tgz is empty or missing: $gzfile"
                        rm -f "$gzfile"
                        failed=$((failed + 1))
                    fi
                else
                    echo "ERROR: Failed to create .tgz for $fitsfile"
                    rm -f "$gzfile"
                    failed=$((failed + 1))
                fi
            else
                echo "Gzipping in-place: $fitsfile -> $gzfile"
                if pigz -p ${SLURM_CPUS_PER_TASK} "$fitsfile"; then
                    processed=$((processed + 1))
                else
                    echo "ERROR: Failed to gzip $fitsfile"
                    failed=$((failed + 1))
                fi
            fi
        else
            echo "ERROR: Failed to backup $fitsfile"
            failed=$((failed + 1))
        fi
    done < "$BATCH_FILE"
else
    echo "ERROR: Unknown batch type: $BATCH_TYPE"
    exit 1
fi

echo ""
echo "Completed at: $(date)"
echo "Processed: $processed | Failed: $failed"

exit $failed
EOFSLURM

chmod +x "$SLURM_SCRIPT"

# ============================================================================
# Submit SLURM jobs
# ============================================================================

echo ""
echo "=== Batch files created in: $WORK_DIR ==="

if $DRY_RUN; then
    echo ""
    echo "=== DRY RUN: Would submit the following jobs ==="
    echo ""
    
    if (( GROUP_BATCHES > 0 )); then
        echo "Group-level batches: $GROUP_BATCHES"
        echo "Command: sbatch --array=0-$((GROUP_BATCHES-1)) $SLURM_SCRIPT group <batch_file> $BACKUP_DIR"
    fi
    
    if (( MEMBER_BATCHES > 0 )); then
        echo "Member-level batches: $MEMBER_BATCHES"
        echo "Command: sbatch --array=0-$((MEMBER_BATCHES-1)) $SLURM_SCRIPT member <batch_file> $BACKUP_DIR"
    fi
    
    echo ""
    echo "To execute, run: bash $0 --execute"
else
    echo ""
    echo "=== Submitting SLURM jobs ==="
    
    # Create submission scripts for each batch type
    if (( GROUP_BATCHES > 0 )); then
        GROUP_SUBMIT="${WORK_DIR}/submit_group.sh"
        cat > "$GROUP_SUBMIT" << EOFSUBMIT
#!/bin/bash
set -euo pipefail

BATCH_LIST="${GROUP_BATCH_LIST}"
BATCH_NUM=\${SLURM_ARRAY_TASK_ID}
BATCH_FILE=\$(sed -n "\$((BATCH_NUM + 1))p" "\$BATCH_LIST")

${SLURM_SCRIPT} group "\$BATCH_FILE" "${BACKUP_DIR}"
EOFSUBMIT
        chmod +x "$GROUP_SUBMIT"
        
        GROUP_JOB=$(sbatch --array=0-$((GROUP_BATCHES-1)) --parsable \
            --job-name=pigz_group \
            --output="${WORK_DIR}/group_%A_%a.log" \
            --time=24:00:00 \
            --ntasks=1 \
            --cpus-per-task=16 \
            --mem=32gb \
            --qos=astronomy-dept-b \
            --account=astronomy-dept \
            "$GROUP_SUBMIT")
        echo "Submitted group-level job array: $GROUP_JOB (${GROUP_BATCHES} tasks)"
    fi
    
    if (( MEMBER_BATCHES > 0 )); then
        MEMBER_SUBMIT="${WORK_DIR}/submit_member.sh"
        cat > "$MEMBER_SUBMIT" << EOFSUBMIT
#!/bin/bash
set -euo pipefail

BATCH_LIST="${MEMBER_BATCH_LIST}"
BATCH_NUM=\${SLURM_ARRAY_TASK_ID}
BATCH_FILE=\$(sed -n "\$((BATCH_NUM + 1))p" "\$BATCH_LIST")

mkdir -p "${BACKUP_DIR}"

${SLURM_SCRIPT} member "\$BATCH_FILE" "${BACKUP_DIR}"
EOFSUBMIT
        chmod +x "$MEMBER_SUBMIT"
        
        MEMBER_JOB=$(sbatch --array=0-$((MEMBER_BATCHES-1)) --parsable \
            --job-name=pigz_member \
            --output="${WORK_DIR}/member_%A_%a.log" \
            --time=24:00:00 \
            --ntasks=1 \
            --cpus-per-task=16 \
            --mem=32gb \
            --qos=astronomy-dept-b \
            --account=astronomy-dept \
            "$MEMBER_SUBMIT")
        echo "Submitted member-level job array: $MEMBER_JOB (${MEMBER_BATCHES} tasks)"
    fi
    
    echo ""
    echo "Monitor jobs with: squeue -u \$USER"
    echo "View logs in: $WORK_DIR"
fi

echo ""
echo "=== Summary ==="
echo "Work directory: $WORK_DIR"
echo "Group batches: $GROUP_BATCHES (${total_files} files)"
echo "Member batches: $MEMBER_BATCHES (${tt_total} files)"