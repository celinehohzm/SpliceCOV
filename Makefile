# Makefile for SpliceCOV (BSD+GNU make compatible, ASCII only)

SHELL   := /usr/bin/env bash
PREFIX  ?= /usr/local
BINDIR  := $(PREFIX)/bin
SHAREDIR:= $(PREFIX)/share/splicecov
PKGNAME := splicecov

# Derive version (fallback if git not available)
VERSION := $(shell git describe --tags --always --dirty 2>/dev/null || echo 0.0.0)

# A default guess; we’ll normalize at install time.
ENTRY   := scripts/splicecov.sh

# Directories to ship into SHAREDIR
SHIP_DIRS := scripts bin

# Optional Python requirements file
REQ     := requirements.txt

.PHONY: release install uninstall check-deps python-deps print-locations help

release: install
	@echo "Installed $(PKGNAME) $(VERSION) to $(PREFIX)"
	@echo "  Launcher:     $(BINDIR)/$(PKGNAME)"
	@echo "  Shared files: $(SHAREDIR)"
	@echo "Run: $(PKGNAME) -h"

install: check-deps python-deps print-locations
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

	@echo "→ Installed files under $(SHAREDIR):"
	@find "$(SHAREDIR)" -maxdepth 2 -type f -print | sed 's/^/   /'

	@echo "Installing launcher to $(BINDIR)/$(PKGNAME)"
	@{ \
	  # Figure out which entry script exists in the source tree
	  printf 'ENTRY_SRC='; \
	  if [ -f scripts/splicecov.sh ]; then \
	    printf '%s\n' 'scripts/splicecov.sh'; \
	  elif [ -f scripts/spliceCOV.sh ]; then \
	    printf '%s\n' 'scripts/spliceCOV.sh'; \
	  else \
	    printf '%s\n' 'scripts/splicecov.sh'; \
	  fi; \
	} > .entry.tmp
	@ENTRY_SRC="$$(cat .entry.tmp)"; rm -f .entry.tmp; \
	{ \
	  printf '%s\n' '#!/usr/bin/env bash'; \
	  printf '%s\n' 'set -euo pipefail'; \
	  printf 'SHAREDIR=%s\n' '$(SHAREDIR)'; \
	  printf 'ENTRY=%s\n' "$$ENTRY_SRC"; \
	  printf 'BIN_LOCAL_ENTRY=%s\n' '$(BINDIR)/.splicecov_entry.sh'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# Export helpers dir for Python scripts (Option B path resolution)'; \
	  printf '%s\n' 'export SPLICECOV_HELPERS_DIR="${SPLICECOV_HELPERS_DIR:-$$SHAREDIR/scripts}"'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# Normalize ENTRY if baked name does not exist (case swap support)'; \
	  printf '%s\n' 'if [[ ! -x "$$SHAREDIR/$$ENTRY" ]]; then'; \
	  printf '%s\n' '  for alt in "scripts/splicecov.sh" "scripts/spliceCOV.sh"; do'; \
	  printf '%s\n' '    if [[ -x "$$SHAREDIR/$$alt" ]]; then ENTRY="$$alt"; break; fi'; \
	  printf '%s\n' '  done'; \
	  printf '%s\n' 'fi'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 1) Try baked/normalized share path'; \
	  printf '%s\n' 'if [[ -x "$$SHAREDIR/$$ENTRY" ]]; then'; \
	  printf '%s\n' '  exec "$$SHAREDIR/$$ENTRY" "$$@"'; \
	  printf '%s\n' 'fi'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 2) Fallbacks relative to launcher'; \
	  printf '%s\n' 'LAUNCHER_DIR="$$(cd "$$(dirname "$${BASH_SOURCE[0]}")" && pwd)"'; \
	  printf '%s\n' 'CANDIDATES=('; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../share/splicecov/$$ENTRY"'; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../../share/splicecov/$$ENTRY"'; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../share/splicecov/scripts/splicecov.sh"'; \
	  printf '%s\n' '  "$$LAUNCHER_DIR/../share/splicecov/scripts/spliceCOV.sh"'; \
	  printf '%s\n' ')'; \
	  printf '%s\n' 'for cand in "$${CANDIDATES[@]}"; do'; \
	  printf '%s\n' '  if [[ -x "$$cand" ]]; then'; \
	  printf '%s\n' '    export SPLICECOV_HELPERS_DIR="$$(cd "$$(dirname "$$cand")" && pwd)";'; \
	  printf '%s\n' '    exec "$$cand" "$$@"'; \
	  printf '%s\n' '  fi'; \
	  printf '%s\n' 'done'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 3) Dev checkout fallback'; \
	  printf '%s\n' 'SCRIPT_DIR_GUESS="$$(cd "$$LAUNCHER_DIR/../scripts" 2>/dev/null && pwd || true)"'; \
	  printf '%s\n' 'for alt in "splicecov.sh" "spliceCOV.sh"; do'; \
	  printf '%s\n' '  if [[ -n "$$SCRIPT_DIR_GUESS" && -x "$$SCRIPT_DIR_GUESS/$$alt" ]]; then'; \
	  printf '%s\n' '    export SPLICECOV_HELPERS_DIR="$$SCRIPT_DIR_GUESS"'; \
	  printf '%s\n' '    exec "$$SCRIPT_DIR_GUESS/$$alt" "$$@"'; \
	  printf '%s\n' '  fi'; \
	  printf '%s\n' 'done'; \
	  printf '%s\n' ''; \
	  printf '%s\n' '# 4) Last resort: local copy beside the launcher'; \
	  printf '%s\n' 'if [[ -x "$$BIN_LOCAL_ENTRY" ]]; then'; \
	  printf '%s\n' '  export SPLICECOV_HELPERS_DIR="$${SPLICECOV_HELPERS_DIR:-$$SHAREDIR/scripts}"'; \
	  printf '%s\n' '  exec "$$BIN_LOCAL_ENTRY" "$$@"'; \
	  printf '%s\n' 'fi'; \
	  printf '%s\n' ''; \
	  printf '%s\n' 'echo "Error: cannot find SpliceCOV entrypoint ($$ENTRY)." >&2'; \
	  printf '%s\n' 'echo "Tried:" >&2'; \
	  printf '%s\n' 'echo "  - $$SHAREDIR/$$ENTRY" >&2'; \
	  printf '%s\n' 'for cand in "$${CANDIDATES[@]}"; do echo "  - $$cand" >&2; done'; \
	  printf '%s\n' 'echo "  - $$BIN_LOCAL_ENTRY" >&2'; \
	  printf '%s\n' 'exit 1'; \
	} > "$(BINDIR)/$(PKGNAME)"
	@chmod +x "$(BINDIR)/$(PKGNAME)"

	@echo "Installing local fallback entry: $(BINDIR)/.splicecov_entry.sh"
	@{ \
	  if [ -f scripts/splicecov.sh ]; then \
	    install -Dm755 "scripts/splicecov.sh" "$(BINDIR)/.splicecov_entry.sh"; \
	  elif [ -f scripts/spliceCOV.sh ]; then \
	    install -Dm755 "scripts/spliceCOV.sh" "$(BINDIR)/.splicecov_entry.sh"; \
	  else \
	    echo "ERROR: neither scripts/splicecov.sh nor scripts/spliceCOV.sh exists." >&2; \
	    exit 1; \
	  fi; \
	}

uninstall:
	@echo "Removing launcher: $(BINDIR)/$(PKGNAME)"
	@rm -f "$(BINDIR)/$(PKGNAME)" || true
	@echo "Removing local fallback entry: $(BINDIR)/.splicecov_entry.sh"
	@rm -f "$(BINDIR)/.splicecov_entry.sh" || true
	@echo "Removing shared files: $(SHAREDIR)"
	@rm -rf "$(SHAREDIR)" || true
	@echo "Uninstalled $(PKGNAME)"

check-deps:
	@echo "Running preflight checks"
	@command -v bash >/dev/null 2>&1   || { echo "Missing: bash" >&2; exit 1; }
	@command -v awk  >/dev/null 2>&1   || { echo "Missing: awk"  >&2; exit 1; }
	@command -v sort >/dev/null 2>&1   || { echo "Missing: sort" >&2; exit 1; }
	@command -v comm >/dev/null 2>&1   || { echo "Missing: comm (coreutils)" >&2; exit 1; }
	@command -v python3 >/dev/null 2>&1|| { echo "Missing: python3" >&2; exit 1; }
	@command -v bigWigToBedGraph >/dev/null 2>&1 || { echo "Missing: bigWigToBedGraph (UCSC utils)" >&2; exit 1; }
	@echo "Core dependencies found"
	@python3 -c "import lightgbm" >/dev/null 2>&1 && echo "Python: lightgbm available" || echo "Python: lightgbm not found (will install if $(REQ) exists)"

python-deps:
	@if [ -f "$(REQ)" ]; then \
	  echo "Installing Python deps from $(REQ) (user scope)"; \
	  python3 -m pip install --user -r "$(REQ)"; \
	else \
	  echo "No $(REQ); skipping Python deps"; \
	fi

print-locations:
	@echo "PKGNAME  = $(PKGNAME)"
	@echo "VERSION  = $(VERSION)"
	@echo "PREFIX   = $(PREFIX)"
	@echo "BINDIR   = $(BINDIR)"
	@echo "SHAREDIR = $(SHAREDIR)"
	@echo "ENTRY    = $(ENTRY)"

help:
	@echo "SpliceCOV Make targets"
	@echo "  make release                         Install into $$PREFIX (default: /usr/local)"
	@echo "  make PREFIX=$$HOME/.local release    Install into user prefix"
	@echo "  make uninstall                       Remove installed launcher and shared dir"
	@echo "  make help                            Show this help"
