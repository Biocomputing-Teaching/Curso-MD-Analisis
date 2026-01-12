OUTDIR=../docs/pdf
DRIVER=episode-driver.tex
SYNC_TALKS=sync_talks_figures.py
EPISODE_SLUGS=01-getting-started 02-running-simulations 03-model-building 04-advanced-examples 05-trajectory-analysis 06-pyemma 07-deeptime
EPISODE_PDFS=$(patsubst %, $(OUTDIR)/curso-md-analisis-ep%.pdf, $(EPISODE_SLUGS))

.PHONY: all episode-pdfs sync-figures clean

all: episode-pdfs sync-figures

episode-pdfs: $(EPISODE_PDFS)

$(OUTDIR)/curso-md-analisis-ep%.pdf: $(DRIVER) course_common.tex or_common.tex master.bib episodes/%.tex FCTE.png
	mkdir -p $(OUTDIR)
	@JOBNAME=$(basename $(notdir $@)) && \
	EPFILE=episodes/$* && \
	SUBTITLE=$$(python ../scripts/episode_subtitle.py $*) && \
	pdflatex -interaction=nonstopmode -halt-on-error -jobname=$$JOBNAME "\def\EpisodeFile{$$EPFILE}\def\EpisodeSubtitle{$$SUBTITLE}\input{$(DRIVER)}" && \
	biber $$JOBNAME && \
	pdflatex -interaction=nonstopmode -halt-on-error -jobname=$$JOBNAME "\def\EpisodeFile{$$EPFILE}\def\EpisodeSubtitle{$$SUBTITLE}\input{$(DRIVER)}" && \
	mv -f $$JOBNAME.pdf $@

sync-figures: $(SYNC_TALKS) curso-md-analisis.tex
	python $(SYNC_TALKS)

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc \
		$(addsuffix .aux,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .log,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .bcf,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .blg,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .bbl,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .run.xml,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(addsuffix .synctex.gz,$(patsubst %,curso-md-analisis-ep%,$(EPISODE_SLUGS))) \
		$(OUTDIR)/curso-md-analisis-ep*.pdf
