.PHONY: all sync-code sync-notebooks standardize-notebooks update-episode-toc render-notebooks check-links check check-sync

# check-links ha de ser el darrer, perquè els sync-* actualitzen
# fitxers que comprova check-links
all:  sync-code sync-notebooks update-episode-toc render-notebooks check-links

sync-code:
	python scripts/sync_episode_code_blocks.py

sync-notebooks:
	python scripts/sync_notebooks_from_scripts.py
	python scripts/standardize_notebook_structure.py

standardize-notebooks:
	python scripts/standardize_notebook_structure.py

update-episode-toc:
	python scripts/update_episode_toc.py

render-notebooks:
	python scripts/render_notebooks_html.py

check-links:
	python scripts/check_links.py

check-sync:
	python scripts/sync_notebooks_from_scripts.py --check

check: sync-code sync-notebooks render-notebooks check-links check-sync
