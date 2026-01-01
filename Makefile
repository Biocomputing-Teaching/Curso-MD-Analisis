.PHONY: all sync-code sync-notebooks

all: sync-notebooks sync-code

sync-code:
	python scripts/sync_episode_code_blocks.py

sync-notebooks:
	python scripts/sync_notebooks_from_scripts.py
