.PHONY: all sync-code sync-notebooks check-links check

# check-links ha de ser el darrer, perquè els sync-* actualitzen
# fitxers que comprova check-links
all:  sync-code sync-notebooks check-links

sync-code:
	python scripts/sync_episode_code_blocks.py

sync-notebooks:
	python scripts/sync_notebooks_from_scripts.py

check-links:
	python scripts/check_links.py
