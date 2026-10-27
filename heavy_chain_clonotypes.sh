#!/usr/bin/env bash
#run in docker
#pwd: file with fasta sequences from paired AbNGS
#docker pull immcantation/suite:4.3.0
#docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.3.0 bash 
# change PREFIX to sequences_heavy, sequences_kappa or sequences_lambda
set -euo pipefail

export IGDATA=/usr/local/share/igblast
igblast_dir=/usr/local/share/igblast

OUT_DIR=/data

PREFIX="sequences_heavy" #change to sequences_kappa or sequences_lambda

mkdir -p igblast_chunks db_chunks "$OUT_DIR"

# ============================================================
# SETTINGS FOR 64-CPU / 500-GB MACHINE
# ============================================================

N_JOBS=8
N_THREADS=4
NPROC=20

# ============================================================
# 1. Split FASTA into ~50k-sequence chunks
# ============================================================
export PREFIX
python3 - <<'PY'
from pathlib import Path
import os

input_file = Path(f"{os.environ['PREFIX']}.fasta")
output_dir = Path("igblast_chunks")
chunk_size = 50_000

output_dir.mkdir(exist_ok=True)

chunk = 0
records = 0
out = None

with input_file.open() as f:
    while True:
        header = f.readline()

        if not header:
            break

        sequence = f.readline()

        if records % chunk_size == 0:
            if out:
                out.close()

            chunk += 1
            path = output_dir / f"chunk_{chunk:04d}.fasta"
            out = path.open("w")

        out.write(header)
        out.write(sequence)

        records += 1

if out:
    out.close()

print(f"Split {records:,} sequences into {chunk} chunks")
PY


# ============================================================
# 2. IgBLAST
# ============================================================

find igblast_chunks -name "*.fasta" -print0 |
xargs -0 -n 1 -P "$N_JOBS" bash -c '
    f="$1"
    base=$(basename "$f" .fasta)

    echo "IgBLAST: $base"

    igblastn \
        -germline_db_V "'"$igblast_dir"'/database/imgt_human_ig_v" \
        -germline_db_D "'"$igblast_dir"'/database/imgt_human_ig_d" \
        -germline_db_J "'"$igblast_dir"'/database/imgt_human_ig_j" \
        -auxiliary_data "'"$igblast_dir"'/optional_file/human_gl.aux" \
        -domain_system imgt \
        -ig_seqtype Ig \
        -organism human \
        -num_threads "'"$N_THREADS"'" \
        -outfmt "7 std qseq sseq btop" \
        -query "$f" \
        -out "igblast_chunks/${base}.fmt7"
' _

# ============================================================
# 3. MakeDb
# ============================================================

find igblast_chunks -name "*.fmt7" -print0 |
xargs -0 -n 1 -P "$N_JOBS" bash -c '
    f="$1"
    base=$(basename "$f" .fmt7)

    echo "MakeDb: $base"

    MakeDb.py igblast \
        -i "$f" \
        -s "igblast_chunks/${base}.fasta" \
        -r /usr/local/share/germlines/imgt/human/vdj/ \
        --extended \
        -o "db_chunks/${base}.tsv"
' _

# ============================================================
# 4. Merge databases
# ============================================================

DB="${PREFIX}_db-pass.tsv"

first=1

for f in db_chunks/*.tsv; do
    if [ "$first" -eq 1 ]; then
        head -n 1 "$f" > "$DB"
        first=0
    fi

    tail -n +2 "$f" >> "$DB"
done

