PY ?= python3
.PHONY: all verify params cover bench clean
all: verify
# Check every number reported in the paper against an independent re-derivation (about 1 minute)
verify:
	$(PY) tests/verify_results.py
# Recompute the stored data from scratch (slow; run from code/ because of relative paths)
params:
	cd code && $(PY) sizes.py
cover:
	cd code && $(PY) cover.py
bench:
	cd code && $(PY) bench.py
clean:
	rm -rf code/__pycache__
