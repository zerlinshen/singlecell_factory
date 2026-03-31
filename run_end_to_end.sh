#!/bin/bash
set -euo pipefail

# End-to-end scRNA-seq Analysis Pipeline
# From Cell Ranger count to Scanpy analysis

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CELLRANGER_DIR="$SCRIPT_DIR/envs/cellranger-10.0.0"
FASTQ_DIR="$SCRIPT_DIR/fastq/pbmc_1k_v3_fastqs"
REFERENCE_DIR="$SCRIPT_DIR/reference/refdata-gex-GRCh38-2024-A"
REPO_OUTPUT_DIR="$SCRIPT_DIR/output"
CONDA_ENV="sc10x"

RUN_ID="pbmc_1k_v3_count"
FASTQ_SAMPLE="pbmc_1k_v3"
EXPECTED_CELLS="1000"
SCANPY_PROJECT_BASE="pbmc_1k_v3_scanpy"
SCANPY_PROJECT=""

SKIP_CELLRANGER=false
SKIP_SCANPY=false
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() { echo -e "${GREEN}[INFO]${NC} $(date +"%H:%M:%S") - $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $(date +"%H:%M:%S") - $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $(date +"%H:%M:%S") - $1"; }

build_timestamped_name() {
    local base_name="$1"
    echo "${base_name}_$(date +%Y%m%d_%H%M%S)"
}

activate_conda() {
    local conda_base
    conda_base="$(conda info --base)"
    # shellcheck disable=SC1090
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
}

check_prerequisites() {
    print_info "Checking prerequisites..."

    if [ ! -x "$CELLRANGER_DIR/bin/cellranger" ]; then
        print_error "Cell Ranger not found at: $CELLRANGER_DIR"
        exit 1
    fi
    if [ ! -f "$REFERENCE_DIR/genes/genes.gtf.gz" ]; then
        print_error "Reference genome not found at: $REFERENCE_DIR"
        exit 1
    fi
    if [ ! -d "$FASTQ_DIR" ]; then
        print_error "FASTQ directory not found: $FASTQ_DIR"
        exit 1
    fi
    if ! conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
        print_error "Conda environment '$CONDA_ENV' not found"
        exit 1
    fi
    mkdir -p "$REPO_OUTPUT_DIR" "$SCRIPT_DIR/.mplconfig" /tmp/numba_cache
    print_info "All prerequisites met"
}

run_cellranger() {
    local cellranger_out="$SCRIPT_DIR/$RUN_ID"
    if [ -d "$cellranger_out/outs/filtered_feature_bc_matrix" ]; then
        print_warning "Cell Ranger output already exists: $cellranger_out"
        print_info "Reusing existing Cell Ranger output"
        return 0
    fi

    print_info "Running Cell Ranger count: $RUN_ID"
    cd "$CELLRANGER_DIR"
    # shellcheck disable=SC1091
    source sourceme.bash
    cellranger count \
        --id="$RUN_ID" \
        --fastqs="$FASTQ_DIR" \
        --sample="$FASTQ_SAMPLE" \
        --transcriptome="$REFERENCE_DIR" \
        --expect-cells="$EXPECTED_CELLS"
    print_info "Cell Ranger completed successfully: $cellranger_out"
}

run_scanpy() {
    local tenx_dir="$SCRIPT_DIR/$RUN_ID/outs/filtered_feature_bc_matrix"
    print_info "Running optimized Scanpy pipeline"

    activate_conda
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    MPLCONFIGDIR="$SCRIPT_DIR/.mplconfig" \
    python "$SCRIPT_DIR/downstream_scanpy/scanpy_pipeline_optimized.py" \
        --tenx_dir="$tenx_dir" \
        --tenx_var_names="gene_symbols" \
        --project="$SCANPY_PROJECT" \
        --base_output_dir="$REPO_OUTPUT_DIR"
    conda deactivate
    print_info "Scanpy pipeline completed successfully"
}

generate_summary() {
    local cellranger_out="$SCRIPT_DIR/$RUN_ID"
    local analysis_out="$REPO_OUTPUT_DIR/$SCANPY_PROJECT"

    print_info "Generating pipeline summary"
    echo "==========================================="
    echo "scRNA-seq Analysis Pipeline Summary"
    echo "==========================================="
    echo "Date: $(date)"
    echo
    echo "Cell Ranger output: $cellranger_out"
    if [ -f "$cellranger_out/outs/web_summary.html" ]; then
        echo "Cell Ranger Web Summary: $cellranger_out/outs/web_summary.html"
    fi
    echo
    echo "Scanpy output: $analysis_out"
    if [ -f "$analysis_out/data/processed.h5ad" ]; then
        echo "Processed Data: $analysis_out/data/processed.h5ad"
    fi
    if [ -f "$analysis_out/tables/cluster_annotation.csv" ]; then
        echo "Cluster Annotation: $analysis_out/tables/cluster_annotation.csv"
    fi
    echo "==========================================="
}

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --help, -h          Show this help message
  --dry-run           Check prerequisites and print planned actions
  --skip-cellranger   Skip Cell Ranger stage
  --skip-scanpy       Skip Scanpy stage
  --project-base=NAME Base Scanpy project name; timestamp is appended automatically
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --dry-run) DRY_RUN=true ;;
        --skip-cellranger) SKIP_CELLRANGER=true ;;
        --skip-scanpy) SKIP_SCANPY=true ;;
        --project-base=*) SCANPY_PROJECT_BASE="${1#*=}" ;;
        *) print_error "Unknown argument: $1"; usage; exit 1 ;;
    esac
    shift
done

main() {
    print_info "Starting scRNA-seq Analysis Pipeline"
    SCANPY_PROJECT="$(build_timestamped_name "$SCANPY_PROJECT_BASE")"
    check_prerequisites

    if [ "$DRY_RUN" = true ]; then
        print_info "Dry run complete. Planned:"
        print_info "  Cell Ranger run id: $RUN_ID"
        print_info "  Scanpy project base: $SCANPY_PROJECT_BASE"
        print_info "  Scanpy project actual: $SCANPY_PROJECT"
        print_info "  Output root: $REPO_OUTPUT_DIR"
        exit 0
    fi

    local CR_OUT="$SCRIPT_DIR/$RUN_ID/outs/filtered_feature_bc_matrix"
    if [ "$SKIP_CELLRANGER" = false ]; then
        run_cellranger
    else
        print_warning "Skipping Cell Ranger stage by user request"
        if [ ! -d "${CR_OUT}" ] || [ ! -f "${CR_OUT}/matrix.mtx.gz" ]; then
            print_error "Cell Ranger output not found at ${CR_OUT}. Cannot skip cellranger without existing data."
            exit 1
        fi
    fi

    if [ "$SKIP_SCANPY" = false ]; then
        run_scanpy
    else
        print_warning "Skipping Scanpy stage by user request"
    fi

    generate_summary
    print_info "Pipeline completed successfully"
}

main
