#!/usr/bin/env bash
set -euo pipefail

########################################
# CONFIG
########################################

eval "$(micromamba shell hook --shell bash)"
micromamba activate bacterial_wgs

export TMPDIR=${TMPDIR:-/tmp}

# Generic paths for seamless Docker volume mounting
export KRAKEN_DB=${KRAKEN_DB:-/databases/kraken2}
export CHECKM_DB=${CHECKM_DB:-/databases/checkm}
export ANTISMASH_DB=${ANTISMASH_DB:-/databases/antismash}

# Container-native paths for CGView
CGVIEW_DIR=/opt/tools/cgview
CGVIEW_JAR=${CGVIEW_JAR:-/opt/tools/cgview.jar}
CGVIEW_XML_BUILDER=${CGVIEW_DIR}/scripts/cgview_xml_builder/cgview_xml_builder.pl

# RESOURCE MANAGEMENT FOR BATCH PROCESSING
MAX_PARALLEL=5      # Processing 5 samples simultaneously
THREADS=6           # 5 samples * 6 threads = 30 total active threads (Safe for node)
INPUT_DIR=$(pwd)

########################################
# TOOL CHECK
########################################
for tool in fastqc fastp spades.py prokka kraken2 checkm multiqc perl mlst abricate macsyfinder parallel blastp blastn; do
    command -v $tool >/dev/null || { echo "ERROR: $tool missing"; exit 1; }
done

# Java detection (base env fallback)
if command -v java >/dev/null; then
    JAVA_CMD="java"
elif micromamba run -n base java -version >/dev/null 2>&1; then
    JAVA_CMD="micromamba run -n base java"
else
    echo "ERROR: java not found"
    exit 1
fi

# CGView Perl module check
if ! micromamba run -n bacterial_wgs perl -MTie::IxHash -e '1' >/dev/null 2>&1; then
    echo "ERROR: Perl module Tie::IxHash missing in bacterial_wgs env."
    exit 1
fi

########################################
# SAMPLE FUNCTION
########################################
run_sample() {

    # Re-initialize micromamba for the isolated subshell
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate bacterial_wgs

    export PATH=$CONDA_PREFIX/bin:$PATH
    export KRAKEN_DB CHECKM_DB ANTISMASH_DB
    export LC_ALL=C LANG=C BLAST_USAGE_REPORT=false

    R1=$1
    R2=$2
    SAMPLE=$(basename "$R1" | sed 's/_R1.*//')

    OUT="${INPUT_DIR}/${SAMPLE}_WGS"

    # Overwrite existing output
    if [[ -d "$OUT" ]]; then
        echo "[${SAMPLE}] Overwriting existing output: removing ${OUT}"
        rm -rf "$OUT"
    fi

    mkdir -p "${OUT}"/{QC/RAW,QC/TRIMMED,TRIM,ASSEMBLY,CHECKM,CHECKM_bins,PROKKA,KRAKEN,MLST,AMR,PLASMID,TXSSCAN,ANTISMASH,PLOTS,TABLES,LOGS,TMP}

    ########################################
    # FASTQC RAW
    ########################################
    echo "[${SAMPLE}] FASTQC (RAW)"
    fastqc -t 2 -o "${OUT}/QC/RAW" "$R1" "$R2" > "${OUT}/LOGS/fastqc_raw.log" 2>&1

    ########################################
    # FASTP & SPADES
    ########################################
    echo "[${SAMPLE}] FASTP"
    fastp -i "$R1" -I "$R2" -o "${OUT}/TRIM/R1.fq.gz" -O "${OUT}/TRIM/R2.fq.gz" --thread ${THREADS} --html "${OUT}/TRIM/fastp.html" > "${OUT}/LOGS/fastp.log" 2>&1

    ########################################
    # FASTQC TRIMMED
    ########################################
    echo "[${SAMPLE}] FASTQC (TRIMMED)"
    fastqc -t 2 -o "${OUT}/QC/TRIMMED" "${OUT}/TRIM/R1.fq.gz" "${OUT}/TRIM/R2.fq.gz" > "${OUT}/LOGS/fastqc_trimmed.log" 2>&1

    echo "[${SAMPLE}] SPADES"
    spades.py --isolate -1 "${OUT}/TRIM/R1.fq.gz" -2 "${OUT}/TRIM/R2.fq.gz" -o "${OUT}/ASSEMBLY" -t ${THREADS} --tmp-dir "${OUT}/TMP" > "${OUT}/LOGS/spades.log" 2>&1

    # Truncate contig names to avoid GenBank LOCUS length errors
    awk '/^>/ {if(seqlen>=1000) print header"\n"seq; header=$0; seq=""; seqlen=0; next} {seq=seq$0; seqlen+=length($0)} END {if(seqlen>=1000) print header"\n"seq}' \
    "${OUT}/ASSEMBLY/contigs.fasta" | sed 's/_length.*//' > "${OUT}/ASSEMBLY/filtered_contigs.fasta"

    CONTIGS="${OUT}/ASSEMBLY/filtered_contigs.fasta"

    ########################################
    # CHECKM
    ########################################
    echo "[${SAMPLE}] CHECKM"
    checkm data setRoot "$CHECKM_DB" > "${OUT}/LOGS/checkm_db.log" 2>&1
    cp "$CONTIGS" "${OUT}/CHECKM_bins/${SAMPLE}.fa"
    checkm lineage_wf "${OUT}/CHECKM_bins" "${OUT}/CHECKM" -x fa -t ${THREADS} --reduced_tree > "${OUT}/LOGS/checkm.log" 2>&1 || true

    if [[ -f "${OUT}/CHECKM/lineage.ms" ]]; then
        checkm qa "${OUT}/CHECKM/lineage.ms" "${OUT}/CHECKM" -o 2 > "${OUT}/TABLES/checkm.tsv" 2>&1 || true
    fi

    ########################################
    # PROKKA (SHADOW BIN BYPASS)
    ########################################
    echo "[${SAMPLE}] PROKKA"
    export TMPDIR="${OUT}/TMP"

    # THE FIX: Create a "Shadow Bin" that contains everything EXCEPT parallel.
    SHADOW_BIN="${OUT}/TMP/shadow_bin"
    mkdir -p "$SHADOW_BIN"
    
    # Symlink all conda binaries into the shadow bin, skipping parallel
    ls "$CONDA_PREFIX/bin" | grep -v '^parallel$' | while read -r bin_name; do
        ln -s "$CONDA_PREFIX/bin/$bin_name" "$SHADOW_BIN/$bin_name" 2>/dev/null || true
    done

    # Strip the original conda bin from PATH to ensure Prokka cannot find parallel anywhere
    CLEAN_PATH=$(echo "$PATH" | sed "s|$CONDA_PREFIX/bin||g" | sed 's/::/:/g' | sed 's/^://' | sed 's/:$//')
    
    # Run Prokka inside this isolated sandbox.
    # It will complain "Could not find 'parallel'", which is exactly what we want.
    PATH="$SHADOW_BIN:$CLEAN_PATH" prokka "$CONTIGS" \
        --outdir "${OUT}/PROKKA" \
        --prefix "${SAMPLE}" \
        --cpus 1 \
        --force \
        > "${OUT}/LOGS/prokka.log" 2>&1 || true

    ########################################
    # KRAKEN2
    ########################################
    echo "[${SAMPLE}] KRAKEN2"
    kraken2 --db "$KRAKEN_DB" --threads ${THREADS} --report "${OUT}/KRAKEN/report.txt" --output "${OUT}/KRAKEN/output.txt" --paired "${OUT}/TRIM/R1.fq.gz" "${OUT}/TRIM/R2.fq.gz" > "${OUT}/LOGS/kraken2.log" 2>&1

    ########################################
    # SANKEY PLOTS
    ########################################
    echo "[${SAMPLE}] SANKEY"
    python <<EOF
import plotly.graph_objects as go
import os

report="${OUT}/KRAKEN/report.txt"

def build_sankey(min_reads=0, out_html="sankey.html"):
    if not os.path.exists(report):
        return

    nodes=[]
    links=[]
    stack=[]

    with open(report) as f:
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)<6:
                continue

            reads=int(p[1])
            name=p[5]
            depth=(len(name)-len(name.lstrip()))//2
            name=name.strip()

            while len(stack)>depth:
                stack.pop()

            if reads >= min_reads:
                if stack:
                    links.append((stack[-1],name,reads))
                nodes.append(name)
                stack.append(name)

    if not links:
        return

    nodes=list(dict.fromkeys(nodes))
    idx={n:i for i,n in enumerate(nodes)}

    fig=go.Figure(go.Sankey(
        node=dict(label=nodes),
        link=dict(
            source=[idx[s] for s,t,v in links],
            target=[idx[t] for s,t,v in links],
            value=[v for s,t,v in links]
        )
    ))
    fig.write_html(out_html)

build_sankey(min_reads=0, out_html="${OUT}/PLOTS/sankey_full.html")
build_sankey(min_reads=100, out_html="${OUT}/PLOTS/sankey_filtered.html")
EOF

    ########################################
    # MLST & AMR
    ########################################
    echo "[${SAMPLE}] MLST & AMR"
    mlst "$CONTIGS" > "${OUT}/MLST/mlst.tsv" 2>&1
    
    # Abricate runs normally because the PATH sandbox was ONLY applied to the Prokka line above
    abricate --db resfinder "$CONTIGS" > "${OUT}/AMR/resfinder.tsv" 2>&1
    abricate --db card "$CONTIGS" > "${OUT}/AMR/card.tsv" 2>&1
    abricate --summary "${OUT}/AMR/"*.tsv > "${OUT}/TABLES/amr_summary.tsv" 2>&1
    abricate --db plasmidfinder "$CONTIGS" > "${OUT}/PLASMID/plasmid.tsv" 2>&1

    ########################################
    # TXSSCAN
    ########################################
    PROKKA_PROTEINS="${OUT}/PROKKA/${SAMPLE}.faa"
    if [[ -f "$PROKKA_PROTEINS" ]]; then
        echo "[${SAMPLE}] TXSSCAN"
        macsyfinder --models TXSScan --sequence-db "$PROKKA_PROTEINS" --out-dir "${OUT}/TXSSCAN" --db-type unordered > "${OUT}/LOGS/txsscan.log" 2>&1 || true
    else
        echo "[${SAMPLE}] SKIPPED TXSSCAN: no Prokka .faa found" | tee -a "${OUT}/LOGS/txsscan.log"
    fi

    ########################################
    # ANTISMASH
    ########################################
    ANTISMASH_IN="${OUT}/PROKKA/${SAMPLE}.gbk"
    ANTISMASH_LOG="${OUT}/LOGS/antismash.log"
    if [[ -f "$ANTISMASH_IN" ]]; then
        echo "[${SAMPLE}] ANTISMASH"
        micromamba run -n antismash_env antismash --taxon bacteria --genefinding-tool none --cpus ${THREADS} --databases "${ANTISMASH_DB}" --output-dir "${OUT}/ANTISMASH" "$ANTISMASH_IN" > "$ANTISMASH_LOG" 2>&1 || true
    else
        echo "[${SAMPLE}] SKIPPED ANTISMASH: no Prokka .gbk found" | tee -a "$ANTISMASH_LOG"
    fi

    ########################################
    # CGVIEW
    ########################################
    CGVIEW_LOG="${OUT}/LOGS/cgview_xml_builder.log"
    if [[ -f "$CGVIEW_JAR" && -f "$CGVIEW_XML_BUILDER" && -f "$ANTISMASH_IN" ]]; then
        echo "[${SAMPLE}] CGVIEW"
        micromamba run -n bacterial_wgs perl "$CGVIEW_XML_BUILDER" -sequence "$ANTISMASH_IN" -gc_content T -gc_skew T -size large-v2 -tick_density 0.05 -draw_divider_rings T -custom showBorder=false title="${SAMPLE}" titleFontSize=200 -output "${OUT}/PLOTS/${SAMPLE}_cgview.xml" > "$CGVIEW_LOG" 2>&1 || true
        $JAVA_CMD -Xmx4g -jar "$CGVIEW_JAR" -i "${OUT}/PLOTS/${SAMPLE}_cgview.xml" -o "${OUT}/PLOTS/${SAMPLE}_cgview.png" -f png > "${OUT}/LOGS/cgview_render.log" 2>&1 || true
    else
        {
            echo "[${SAMPLE}] SKIPPED CGVIEW."
        } | tee -a "$CGVIEW_LOG"
    fi

    ########################################
    # MULTIQC
    ########################################
    echo "[${SAMPLE}] MULTIQC"
    multiqc "${OUT}" -o "${OUT}" --force > "${OUT}/LOGS/multiqc.log" 2>&1

    echo "DONE: $SAMPLE"
}

export CGVIEW_DIR CGVIEW_JAR CGVIEW_XML_BUILDER JAVA_CMD
export KRAKEN_DB CHECKM_DB ANTISMASH_DB INPUT_DIR THREADS CONDA_PREFIX
export -f run_sample

########################################
# AUTO DETECT & RUN
########################################
find "${INPUT_DIR}" -name "*_R1*.fastq.gz" | sort | while read R1
do
    R2=$(echo "$R1" | sed 's/_R1/_R2/')
    if [[ -f "$R2" ]]; then
        echo "$R1 $R2"
    fi
done | xargs -n 2 -P ${MAX_PARALLEL} bash -c 'run_sample "$1" "$2"' _