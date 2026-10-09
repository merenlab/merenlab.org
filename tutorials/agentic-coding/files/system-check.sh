#!/usr/bin/env bash
# ==============================================================================
# Pre-Workshop Check Checks
# ==============================================================================

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
BOLD='\033[1m'
NC='\033[0m' # No Color

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

print_header() {
    echo -e "\n${BOLD}====================================================${NC}"
    echo -e "${BOLD} $1 ${NC}"
    echo -e "${BOLD}====================================================${NC}"
}

report_pass() {
    echo -e "  [${GREEN} PASS ${NC}] $1"
    ((PASS_COUNT++))
}

report_warn() {
    echo -e "  [${YELLOW} WARN ${NC}] $1"
    echo -e "           ${YELLOW}--> Note:${NC} $2"
    ((WARN_COUNT++))
}

report_fail() {
    echo -e "  [${RED} FAIL ${NC}] $1"
    echo -e "           ${RED}--> Action required:${NC} $2"
    ((FAIL_COUNT++))
}

# Helper: Detect System RAM in GB
get_total_ram_gb() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        local ram_bytes
        ram_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
        echo $((ram_bytes / 1024 / 1024 / 1024))
    elif [ -f /proc/meminfo ]; then
        local ram_kb
        ram_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        echo $((ram_kb / 1024 / 1024))
    else
        echo 0
    fi
}

# Helper: Detect GPU Acceleration
get_gpu_info() {
    # Apple Silicon (Unified Memory / Metal)
    if [[ "$OSTYPE" == "darwin"* ]] && [[ "$(uname -m)" == "arm64" ]]; then
        echo "Apple Silicon (Metal Unified Memory)"
        return
    fi

    # NVIDIA GPU via nvidia-smi
    if command -v nvidia-smi &>/dev/null; then
        local nv_gpu
        nv_gpu=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -n 1)
        if [ -n "$nv_gpu" ]; then
            echo "NVIDIA GPU ($nv_gpu)"
            return
        fi
    fi

    # AMD ROCm
    if command -v rocm-smi &>/dev/null; then
        echo "AMD GPU (ROCm detected)"
        return
    fi

    echo "None detected (CPU execution only)"
}

print_header "1. Checking Core Command-Line Tools & Git"

if command -v git &>/dev/null; then
    GIT_VER=$(git --version)
    report_pass "git is installed ($GIT_VER)"

    GIT_USER=$(git config --global user.name)
    GIT_EMAIL=$(git config --global user.email)

    if [ -n "$GIT_USER" ] && [ -n "$GIT_EMAIL" ]; then
        report_pass "git identity configured ($GIT_USER <$GIT_EMAIL>)"
    else
        report_fail "git identity is missing" "Run 'git config --global user.name \"Your Name\"' and 'git config --global user.email \"your@email.com\"'"
    fi
else
    report_fail "git is not installed" "Install git using your package manager or from https://git-scm.com/"
fi

if command -v curl &>/dev/null; then
    report_pass "curl is available"
else
    report_fail "curl is not installed" "Install curl using your operating system package manager."
fi

print_header "2. Checking OpenCode (Agent)"

if command -v opencode &>/dev/null; then
    OPENCODE_VER=$(opencode --version 2>&1)
    if [[ "$OPENCODE_VER" == *"1.18.33"* ]]; then
        report_pass "OpenCode is installed with exact workshop version ($OPENCODE_VER)"
    else
        report_warn "OpenCode version mismatch ($OPENCODE_VER)" "Workshop targets v1.18.33. Most features will work, but pinned version is recommended."
    fi
else
    report_fail "OpenCode is not installed" "Run: curl -fsSL https://opencode.ai/install | bash -s -- --version 1.18.33"
fi

print_header "3. Checking R & Tidyverse Environment"

if command -v Rscript &>/dev/null; then
    R_VER=$(Rscript -e 'cat(R.version.string)')
    report_pass "R is installed ($R_VER)"

    if Rscript -e 'suppressPackageStartupMessages(library(tidyverse))' &>/dev/null; then
        report_pass "tidyverse R package is loadable"
    else
        report_fail "tidyverse package not found in R" "Open R and run: install.packages(\"tidyverse\")"
    fi
else
    report_fail "R executable (Rscript) not found" "Install R from CRAN or activate the appropriate conda environment (such as anvio-dev :))."
fi

print_header "4. Checking Anvi'o Development Environment"

if command -v anvi-self-test &>/dev/null; then
    ANVIO_VER=$(anvi-self-test -v 2>&1 | head -n 1)
    report_pass "anvi'o dev environment active ($ANVIO_VER)"
else
    report_warn "anvi'o development environment not active/installed" "Ensure 'conda activate anvio-dev' works IF you will follow anvi'o examples."
fi

print_header "5. Checking Ollama Service & Local Models"

if command -v ollama &>/dev/null; then
    OLLAMA_VER=$(ollama --version 2>&1)
    report_pass "Ollama binary is installed ($OLLAMA_VER)"

    if ollama list &>/dev/null; then
        report_pass "Ollama daemon is running"

        PULL_TEST=$(ollama list | grep -E "qwen|minicpm|gpt-oss")
        if [ -n "$PULL_TEST" ]; then
            report_pass "Found at least one workshop model pulled in Ollama"
        else
            report_warn "No workshop models pulled yet" "Run 'ollama pull openbmb/minicpm5-2b' to download the first test model."
        fi
    else
        report_warn "Ollama daemon is not responding" "Start Ollama in another terminal window using 'ollama serve' or launch the desktop app."
    fi
else
    report_fail "Ollama is not installed" "Install Ollama from https://ollama.com/download"
fi

print_header "6. Checking Hardware Resources (RAM & GPU Acceleration)"

TOTAL_RAM=$(get_total_ram_gb)
GPU_INFO=$(get_gpu_info)

# Memory Capacity Evaluation
if [ "$TOTAL_RAM" -ge 22 ]; then
    report_pass "System RAM: ${TOTAL_RAM} GB (Sufficient for larger models like gpt-oss:20b)"
elif [ "$TOTAL_RAM" -ge 16 ]; then
    report_warn "System RAM: ${TOTAL_RAM} GB" "Do not go beyond qwen3.5:4b or minicpm5-2b."
elif [ "$TOTAL_RAM" -ge 8 ]; then
    report_warn "System RAM: ${TOTAL_RAM} GB" "Do not use anything other than minicpm5-2b or equivalent."
elif [ "$TOTAL_RAM" -gt 0 ]; then
    report_fail "System RAM: ${TOTAL_RAM} GB" "Insufficient memory (< 8 GB) to run local models."
else
    report_warn "Could not determine system RAM" "Ensure you have at least 8 GB RAM to run lightweight local models."
fi

# GPU Acceleration Evaluation
if [[ "$GPU_INFO" != *"None detected"* ]]; then
    report_pass "GPU Acceleration: $GPU_INFO"
else
    report_warn "GPU Acceleration: None detected" "Running local models on CPU will be extremely slow. Do not try running models live—focus on the live content or tutorial demonstrations instead."
fi

print_header "Diagnostic Summary"

echo -e "  Passed: ${GREEN}${PASS_COUNT}${NC}"
echo -e "  Warnings: ${YELLOW}${WARN_COUNT}${NC}"
echo -e "  Failures: ${RED}${FAIL_COUNT}${NC}\n"

if [ $FAIL_COUNT -eq 0 ] && [ $WARN_COUNT -eq 0 ]; then
    echo -e "${GREEN}${BOLD}Everything looks perfect! You can raise your thumbs-up for Meren to see.${NC}\n"
elif [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${YELLOW}${BOLD}Not great, not terrible: Core requirements are met, but you should review the hardware/environment warnings.${NC}\n"
else
    echo -e "${RED}${BOLD}Some critical dependencies are missing :/ Please raise items with [ FAIL ] with Meren.${NC}\n"
fi
