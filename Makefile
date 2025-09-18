# Makefile for SpliceCOV
# Usage:
#   make release                # install to /usr/local by default
#   make PREFIX=/opt release    # custom install prefix
#   make uninstall              # remove installed files
#   make help                   # show targets

.RECIPEPREFIX := >
SHELL   := /usr/bin/env bash
PREFIX  ?= /usr/local
BINDIR  := $(PREFIX)/bin
SHAREDIR:= $(PREFIX)/share/splicecov
PKGNAME := splicecov
VERSION := $(shell git describe --tags --always --dirty 2>/dev/null || echo 0.0.0)
ENTRY   := scripts/splicecov.sh
SHIP_DIRS := scripts bin
REQ     := requirements.txt

.PHONY: release install uninstall check-deps python-deps print-locations help

release: install
> echo "✅ Installed $(PKGNAME) $(VERSION) to $(PREFIX)"
> echo "• Launcher:     $(BINDIR)/$(PKGNAME)"
> echo "• Shared files: $(SHAREDIR)"
> echo "Run: $(PKGNAME) -h"

install: check-deps python-deps print-locations
> echo "→ Creating directories"
> mkdir -p "$(BINDIR)" "$(SHAREDIR)"
> echo "→ Copying project files to $(SHAREDIR)"
> for d in $(SHIP_DIRS); do \
>   if [[ -d $$d ]]; then \
>     rsync -a --delete $$d "$(SHAREDIR)/"; \
>   fi; \
> done
> echo "→ Marking scripts executable"
> find "$(SHAREDIR)/scripts" -type f -name "*.sh" -exec chmod +x {} \; || true
> find "$(SHAREDIR)/bin"     -type f -perm -u=x     -exec chmod +x {} \; || true
> echo "→ Installing launcher to $(BINDIR)/$(PKGNAME)"
> cat > "$(BINDIR)/$(PKGNAME)" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
SHAREDIR="{{SHAREDIR}}"
ENTRY="{{ENTRY}}"
export SPLICECOV_HELPERS_DIR="${SPLICECOV_HELPERS_DIR:-$SHAREDIR/scripts}"
if [[ -x "$SHAREDIR/$ENTRY" ]]; then
  exec "$SHAREDIR/$ENTRY" "$@"
elif [[ -x "$ENTRY" ]]; then
  exec "$ENTRY" "$@"
else
  echo "Error: cannot find SpliceCOV entrypoint ($ENTRY). Reinstall or fix launcher." >&2
  exit 1
fi
EOF
> # macOS-safe in-place edit
> sed -i.bak "s#{{SHAREDIR}}#$(SHAREDIR)#g; s#{{ENTRY}}#$(ENTRY)#g" "$(BINDIR)/$(PKGNAME)" && rm -f "$(BINDIR)/$(PKGNAME).bak"
> chmod +x "$(BINDIR)/$(PKGNAME)"

uninstall:
> echo "→ Removing launcher: $(BINDIR)/$(PKGNAME)"
> rm -f "$(BINDIR)/$(PKKNAME)" || true

# fix typo in uninstall path
uninstall:
> echo "→ Removing launcher: $(BINDIR)/$(PKGNAME)"
> rm -f "$(BINDIR)/$(PKGNAME)" || true
> echo "→ Removing shared files: $(SHAREDIR)"
> rm -rf "$(SHAREDIR)" || true
> echo "✅ Uninstalled $(PKGNAME)"

check-deps:
> echo "→ Running preflight checks"
> hash bash 2>/dev/null || (echo "Missing: bash" >&2; exit 1)
> hash awk  2>/dev/null || (echo "Missing: awk"  >&2; exit 1)
> hash sort 2>/dev/null || (echo "Missing: sort" >&2; exit 1)
> hash comm 2>/dev/null || (echo "Missing: comm (coreutils)" >&2; exit 1)
> hash python3 2>/dev/null || (echo "Missing: python3" >&2; exit 1)
> hash rsync 2>/dev/null || (echo "Missing: rsync" >&2; exit 1)
> hash bigWigToBedGraph 2>/dev/null || (echo "Missing: bigWigToBedGraph (UCSC utils)" >&2; exit 1)
> echo "✓ core deps found"
> if python3 -c "import lightgbm" 2>/dev/null; then \
>   echo "✓ Python: lightgbm available"; \
> else \
>   echo "⚠︎ Python: lightgbm not found (will try to install if $(REQ) exists)"; \
> fi

python-deps:
> if [[ -f "$(REQ)" ]]; then \
>   echo "→ Installing Python deps from $(REQ) (user scope)"; \
>   python3 -m pip install --user -r "$(REQ)"; \
> else \
>   echo "→ No $(REQ) found; skipping Python deps"; \
> fi

print-locations:
> echo "PKGNAME  = $(PKGNAME)"
> echo "VERSION  = $(VERSION)"
> echo "PREFIX   = $(PREFIX)"
> echo "BINDIR   = $(BINDIR)"
> echo "SHAREDIR = $(SHAREDIR)"
> echo "ENTRY    = $(ENTRY)"

help:
> echo "SpliceCOV Make targets"
> echo "  make release              Install into \$$PREFIX (default: /usr/local)"
> echo "  make PREFIX=/opt release  Install into custom prefix"
> echo "  make uninstall            Remove installed launcher and shared dir"
> echo "  make help                 Show this help"
