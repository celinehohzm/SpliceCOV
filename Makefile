# Makefile for SpliceCOV (BSD+GNU make compatible, ASCII only)

SHELL     := /usr/bin/env bash

# -------- Paths / names --------
PREFIX    ?= /usr/local
BINDIR    := $(PREFIX)/bin
SHAREDIR  := $(PREFIX)/share/splicecov
PKGNAME   := splicecov

# Python interpreter & minimum version
PYTHON        ?= python3
MIN_PY_MAJOR  := 3
MIN_PY_MINOR  := 10

# Entrypoint inside repo
ENTRY     := scripts/spliceCOV.sh
SHIP_DIRS := scripts bin
REQ       := requirements.txt

# -------- C build (process_tiebrush) --------
CC      ?= cc
CFLAGS  ?= -O2 -std=c17 -Wall -Wextra
LDFLAGS ?=
LDLIBS  ?= -lm

C_SRC   := scripts/process_tiebrush.c
C_BIN   := bin/process_tiebrush

# -------- Version (fallback if git not available) --------
VERSION   := $(shell git describe --tags --always --dirty 2>/dev/null || echo 0.0.0)

# -------- Auto-deps controls --------
SKIP_AUTO_DEPS ?= 0          # set to 1 to disable auto-install of UCSC tools
FORCE_REINSTALL ?= 0         # set to 1 to force re-install of UCSC tool
UCSC_TOOL ?= bigWigToBedGraph

.PHONY: release install uninstall check-deps python-deps build-c clean print-locations help \
        install-ucsc-bw2bg _copy-tree _make-launcher

# =========================================================
# Top-level targets
# =========================================================
release: install
	@echo "Installed $(PKGNAME) $(VERSION) to $(PREFIX)"
	@echo "  Launcher:     $(BINDIR)/$(PKGNAME)"
	@echo "  Shared files: $(SHAREDIR)"
	@echo "Run: $(PKGNAME) -h"

# NOTE: build-c compiles scripts/process_tiebrush.c -> bin/process_tiebrush
# Then _copy-tree ships bin/ into $(SHAREDIR)/bin because SHIP_DIRS includes "bin".
install: check-deps python-deps build-c print-locations _copy-tree _make-launcher
	@echo "→ Installed files under $(SHAREDIR):"
	@find "$(SHAREDIR)" -maxdepth 2 -type f -print | sed 's/^/   /'

uninstall:
	@echo "Removing launcher: $(BINDIR)/$(PKGNAME)"
	@rm -f "$(BINDIR)/$(PKGNAME)" || true
	@echo "Removing shared files: $(SHAREDIR)"
	@rm -rf "$(SHAREDIR)" || true
	@echo "Uninstalled $(PKGNAME)"

# =========================================================
# Dependency checks (with auto-install for UCSC tool)
# =========================================================
check-deps:
	@echo "Running preflight checks"
	@command -v bash    >/dev/null 2>&1 || { echo "Missing: bash"  >&2; exit 1; }
	@command -v awk     >/dev/null 2>&1 || { echo "Missing: awk"   >&2; exit 1; }
	@command -v sort    >/dev/null 2>&1 || { echo "Missing: sort"  >&2; exit 1; }
	@command -v comm    >/dev/null 2>&1 || { echo "Missing: comm (coreutils)" >&2; exit 1; }
	@command -v $(PYTHON) >/dev/null 2>&1 || { echo "Missing: $(PYTHON)" >&2; exit 1; }
	@command -v $(CC)   >/dev/null 2>&1 || { echo "Missing: $(CC) (C compiler)" >&2; exit 1; }
	@# Enforce Python >= 3.10 (portable one-liner; no heredoc)
	@$(PYTHON) -c 'import sys; req=($(MIN_PY_MAJOR),$(MIN_PY_MINOR)); cur=sys.version_info; \
ok=(cur.major,cur.minor)>=req; \
sys.exit(0) if ok else (sys.stderr.write(f"ERROR: Python >= {req[0]}.{req[1]} required; found {cur.major}.{cur.minor}.{cur.micro}\nTip:\n  conda create -n splicecov -c conda-forge -c bioconda python=3.11 lightgbm ucsc-bigwigtobedgraph -y\n  conda activate splicecov\n"), sys.exit(2))'
	@# Try to ensure bigWigToBedGraph exists; attempt auto-install unless opted out
	@if ! command -v $(UCSC_TOOL) >/dev/null 2>&1; then \
	  if [ "$(SKIP_AUTO_DEPS)" != "1" ]; then \
	    echo "Missing: $(UCSC_TOOL) (UCSC utils) — attempting auto-install..."; \
	    $(MAKE) install-ucsc-bw2bg || true; \
	  fi; \
	fi
	@if ! command -v $(UCSC_TOOL) >/dev/null 2>&1; then \
	  echo "ERROR: $(UCSC_TOOL) still missing."; \
	  echo "       Try: conda install -y -c bioconda ucsc-bigwigtobedgraph"; \
	  echo "       Or ensure it is on PATH or set PREFIX/BINDIR appropriately."; \
	  exit 1; \
	fi
	@echo "Core dependencies found"
	@$(PYTHON) -c "import lightgbm" >/dev/null 2>&1 && echo "Python: lightgbm available" || echo "Python: lightgbm not found (will install if $(REQ) exists)"

# Install UCSC bigWigToBedGraph (prefers conda; falls back to UCSC prebuilt)
install-ucsc-bw2bg:
	@set -euo pipefail; \
	if command -v $(UCSC_TOOL) >/dev/null 2>&1 && [ "$(FORCE_REINSTALL)" != "1" ]; then \
	  echo "$(UCSC_TOOL) already present; skip (set FORCE_REINSTALL=1 to force)."; \
	  exit 0; \
	fi; \
	mkdir -p "$(BINDIR)"; \
	if command -v conda >/dev/null 2>&1; then \
	  echo "Trying conda (bioconda) for $(UCSC_TOOL)..."; \
	  if conda install -y -c bioconda ucsc-bigwigtobedgraph >/dev/null 2>&1; then \
	    if command -v $(UCSC_TOOL) >/dev/null 2>&1; then \
	      echo "Installed via conda."; \
	      exit 0; \
	    fi; \
	  fi; \
	  echo "Conda install did not place $(UCSC_TOOL) on PATH, trying UCSC binary download..."; \
	fi; \
	uname_s="$$(uname -s)"; uname_m="$$(uname -m)"; \
	case "$$uname_s:$$uname_m" in \
	  Linux:x86_64)  ucsc_os="linux.x86_64" ;; \
	  Darwin:arm64)  ucsc_os="macOSX.arm64" ;; \
	  Darwin:x86_64) ucsc_os="macOSX.x86_64" ;; \
	  *) echo "Unsupported platform: $$uname_s $$uname_m" >&2; exit 1 ;; \
	esac; \
	url="https://hgdownload.soe.ucsc.edu/admin/exe/$$ucsc_os/$(UCSC_TOOL)"; \
	echo "Downloading $(UCSC_TOOL) from $$url -> $(BINDIR)/$(UCSC_TOOL)"; \
	tmp="$(BINDIR)/$(UCSC_TOOL).tmp.$$RANDOM"; \
	curl -fL --retry 3 -o "$$tmp" "$$url"; \
	chmod +x "$$tmp"; \
	mv "$$tmp" "$(BINDIR)/$(UCSC_TOOL)"; \
	if ! "$(BINDIR)/$(UCSC_TOOL)" 2>&1 | head -n1 >/dev/null; then \
	  echo "Downloaded $(UCSC_TOOL) appears not to run. Aborting." >&2; \
	  exit 1; \
	fi; \
	echo "Installed $(UCSC_TOOL) -> $(BINDIR)/$(UCSC_TOOL)"

# =========================================================
# Python deps
# =========================================================
python-deps:
	@if [ -f "$(REQ)" ]; then \
	  echo "Installing Python deps from $(REQ) using $(PYTHON)"; \
	  $(PYTHON) -m pip install --upgrade pip; \
	  if [ -n "$$CONDA_PREFIX" ] || [ -n "$$VIRTUAL_ENV" ]; then \
	    echo "(detected env) installing into the active environment"; \
	    $(PYTHON) -m pip install -r "$(REQ)"; \
	  else \
	    echo "(no env) installing with --user"; \
	    $(PYTHON) -m pip install --user -r "$(REQ)"; \
	  fi; \
	else \
	  echo "No $(REQ); skipping Python deps"; \
	fi

# =========================================================
# C build
# =========================================================
build-c: $(C_BIN)

$(C_BIN): $(C_SRC)
	@echo "Compiling $< -> $@"
	@mkdir -p "$(@D)"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o "$@" "$<" $(LDLIBS)

clean:
	@rm -f "$(C_BIN)"

# =========================================================
# Install payload + launcher
# =========================================================
_copy-tree:
	@echo "Creating directories"
	@mkdir -p "$(BINDIR)" "$(SHAREDIR)"
	@echo "Copying project files to $(SHAREDIR)"
	@for d in $(SHIP_DIRS); do \
	  if [ -d "$$d" ]; then \
	    if command -v rsync >/dev/null 2>&1; then \
	      rsync -a --delete "$$d" "$(SHAREDIR)/"; \
	    else \
	      rm -rf "$(SHAREDIR)/$${d##*/}" 2>/dev/null || true; \
	      cp -R "$$d" "$(SHAREDIR)/"; \
	    fi; \
	  fi; \
	done
	@echo "Marking scripts executable"
	@find "$(SHAREDIR)/scripts" -type f -name "*.sh" -exec chmod +x {} \; 2>/dev/null || true
	@find "$(SHAREDIR)/bin"     -type f -perm -u=x     -exec chmod +x {} \; 2>/dev/null || true

_make-launcher:
	@echo "Installing launcher to $(BINDIR)/$(PKGNAME)"
	@{ \
	  printf '%s\n' '#!/usr/bin/env bash'; \
	  printf '%s\n' 'set -euo pipefail'; \
	  printf 'SHAREDIR=%s\n' '$(SHAREDIR)'; \
	  printf 'ENTRY=%s\n' '$(ENTRY)'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# Export helpers dir for Python scripts (used by Option B path resolution)'; \
	  printf '%s\n' 'export SPLICECOV_HELPERS_DIR="${SPLICECOV_HELPERS_DIR:-$$SHAREDIR/scripts}"'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 1) Try baked-in path first'; \
	  printf '%s\n' 'if [[ -x "$$SHAREDIR/$$ENTRY" ]]; then'; \
	  printf '%s\n' '  exec "$$SHAREDIR/$$ENTRY" "$$@"'; \
	  printf '%s\n' 'fi'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 2) Fallback: derive share dir from launcher location'; \
	  printf '%s\n' 'LAUNCHER_DIR="$$(cd "$$(dirname "$${BASH_SOURCE[0]}")" && pwd)"'; \
	  printf '%s\n' 'CANDIDATES=('; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../share/splicecov/$$ENTRY"'; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../../share/splicecov/$$ENTRY"'; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../share/splicecov/scripts/spliceCOV.sh"'; \
	  printf '%s\n' ')'; \
	  printf '%s\n' 'for cand in "$${CANDIDATES[@]}"; do'; \
	  printf '%s\n' '  if [[ -x "$$cand" ]]; then'; \
	  printf '%s\n' '    export SPLICECOV_HELPERS_DIR="$$(cd "$$(dirname "$$cand")" && pwd)";'; \
	  printf '%s\n' '    exec "$$cand" "$$@"'; \
	  printf '%s\n' '  fi'; \
	  printf '%s\n' 'done'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 3) Last chance: dev checkout relative paths'; \
	  printf '%s\n' 'SCRIPT_DIR_GUESS="$$(cd "$$LAUNCHER_DIR/../scripts" 2>/dev/null && pwd || true)"'; \
	  printf '%s\n' 'if [[ -n "$$SCRIPT_DIR_GUESS" && -x "$$SCRIPT_DIR_GUESS/spliceCOV.sh" ]]; then'; \
	  printf '%s\n' '  export SPLICECOV_HELPERS_DIR="$$SCRIPT_DIR_GUESS"'; \
	  printf '%s\n' '  exec "$$SCRIPT_DIR_GUESS/spliceCOV.sh" "$$@"'; \
	  printf '%s\n' 'fi'; \
	  printf '%s\n' ''; \
	  printf '%s\n' 'echo "Error: cannot find SpliceCOV entrypoint ($$ENTRY)." >&2'; \
	  printf '%s\n' 'echo "Tried:" >&2'; \
	  printf '%s\n' 'echo "  - $$SHAREDIR/$$ENTRY" >&2'; \
	  printf '%s\n' 'for cand in "$${CANDIDATES[@]}"; do echo "  - $$cand" >&2; done'; \
	  printf '%s\n' 'exit 1'; \
	} > "$(BINDIR)/$(PKGNAME)"
	@chmod +x "$(BINDIR)/$(PKGNAME)"

# =========================================================
# Info / help
# =========================================================
print-locations:
	@echo "PKGNAME  = $(PKGNAME)"
	@echo "VERSION  = $(VERSION)"
	@echo "PREFIX   = $(PREFIX)"
	@echo "BINDIR   = $(BINDIR)"
	@echo "SHAREDIR = $(SHAREDIR)"
	@echo "ENTRY    = $(ENTRY)"
	@echo "PYTHON   = $(PYTHON)"
	@echo "MIN_PY   = $(MIN_PY_MAJOR).$(MIN_PY_MINOR)"
	@echo "CC       = $(CC)"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "C_SRC    = $(C_SRC)"
	@echo "C_BIN    = $(C_BIN)"
	@echo "SKIP_AUTO_DEPS = $(SKIP_AUTO_DEPS)"
	@echo "FORCE_REINSTALL = $(FORCE_REINSTALL)"

help:
	@echo "SpliceCOV Make targets"
	@echo "  make release                          Install into \$$PREFIX (default: /usr/local)"
	@echo "  make PREFIX=\$$HOME/.local release     Install into user prefix"
	@echo "  make build-c                          Compile scripts/process_tiebrush.c -> bin/process_tiebrush"
	@echo "  make clean                            Remove compiled C binary"
	@echo "  make uninstall                        Remove installed launcher and shared dir"
	@echo "  make help                             Show this help"
	@echo ""
	@echo "Variables:"
	@echo "  PYTHON=$(PYTHON) (override to choose interpreter)"
	@echo "  MIN_PY=$(MIN_PY_MAJOR).$(MIN_PY_MINOR) (required minimum)"
	@echo "  CC=$(CC) (override compiler)"
	@echo "  CFLAGS=$(CFLAGS) (override C flags)"
	@echo ""
	@echo "Auto-deps:"
	@echo "  SKIP_AUTO_DEPS=1   Disable auto-install of UCSC tool"
	@echo "  FORCE_REINSTALL=1  Force reinstall of UCSC tool"

