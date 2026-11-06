#!/usr/bin/env bash
set -euo pipefail

source /path/to/user_workspace/miniconda3/etc/profile.d/conda.sh
conda activate nf-core-env

export LANG=C.utf8
export LC_ALL=C.utf8
export LC_CTYPE=C.utf8
export LC_MESSAGES=C.utf8

export JAVA_TOOL_OPTIONS="--add-opens=java.base/java.util=ALL-UNNAMED --add-opens=java.base/java.lang=ALL-UNNAMED --add-opens=java.base/java.io=ALL-UNNAMED -Dfile.encoding=UTF-8 -Dsun.jnu.encoding=UTF-8 ${JAVA_TOOL_OPTIONS:-}"

BASE=/path/to/user_workspace/nf-core_sarek
WF=${BASE}/workflow/sarek-3.8.1
CACHE=${BASE}/singularity_cache
NXF_JAR=/path/to/user_workspace/.nextflow/framework/25.10.4/nextflow-25.10.4-one.jar

export NXF_HOME=/path/to/user_workspace/.nextflow
export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
export NXF_ANSI_LOG=false
export NXF_SINGULARITY_CACHEDIR=${CACHE}
export NXF_APPTAINER_CACHEDIR=${CACHE}

export TMPDIR=${BASE}/tmp
export SINGULARITY_CACHEDIR=${BASE}/tmp/singularity_cache
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR" logs

WES_REGION=/path/to/references/wholeExonBed/company/company_use_agilent_region.bed
GERMLINE_RESOURCE=/path/to/user_workspace/reference/gatk_mutect2_resource/b37/af-only-gnomad.raw.sites.b37.vcf.gz
GERMLINE_RESOURCE_TBI=/path/to/user_workspace/reference/gatk_mutect2_resource/b37/af-only-gnomad.raw.sites.b37.vcf.gz.tbi
RUN_DIR=/path/to/user_workspace/data_WES/sarek_wes_somatic_mutect2_strelka_freebayes_vep_from_fastq
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

for f in \
  "$RUN_DIR/clinical_calling.config" \
  "$RUN_DIR/cohort_samplesheet_wes_somatic_main_from_fastq.csv" \
  "$WES_REGION" \
  "$GERMLINE_RESOURCE" \
  "$GERMLINE_RESOURCE_TBI"
do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing or empty file: $f" >&2
    exit 1
  fi
done

java -jar "$NXF_JAR" -version

java -jar "$NXF_JAR" \
  -c "$RUN_DIR/clinical_calling.config" \
  run "$WF" \
  -profile singularity \
  --input "$RUN_DIR/cohort_samplesheet_wes_somatic_main_from_fastq.csv" \
  --outdir "$RUN_DIR/cohort_results_wes_somatic_mutect2_strelka_freebayes_vep_from_fastq" \
  --tools mutect2,muse,strelka,freebayes,vep \
  --wes \
  --germline_resource "$GERMLINE_RESOURCE" \
  --germline_resource_tbi "$GERMLINE_RESOURCE_TBI" \
  --intervals "$WES_REGION" \
  --save_mapped \
  -work-dir "$RUN_DIR/cohort_work_wes_somatic_mutect2_strelka_freebayes_vep_from_fastq" \
  -ansi-log false \
  -resume
