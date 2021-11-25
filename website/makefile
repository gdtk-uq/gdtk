
DGD_DOC ?= ../../dgd-docs
EILMER_DIR := content/docs/eilmer
L1D_DIR := content/docs/l1d
ESTCN_DIR := content/docs/estcn
NENZF1D_DIR := content/docs/nenzf1d
PYTHON_DIR := content/docs/python

REF_MANS := eilmer-reference-manual \
	l1d-reference-manual \
	estcn-manual \
	nenzf1d-manual \
	gas-dynamic-manual \
	geometry-manual \
	nm-manual

USER_GUIDES := eilmer-user-guide \
	geometry-user-guide \
	gas-user-guide \
	reacting-gas-guide \
	foammesh-user-guide

TECH_NOTES := hybrid-fluxes

all: $(REF_MANS) eilmer-catalog eilmer-faq presentations theses $(USER_GUIDES) $(TECH_NOTES)
	hugo

sync:
	rsync --exclude=pdfs -zav public/ gdtk:/var/www/htdocs/

sync-pdfs:
	rsync -zav public/pdfs/ gdtk:/var/www/htdocs/pdfs/

clean:
	- rm -rf $(EILMER_DIR)/eilmer-catalog-for-hugo/gallery
	- rm -rf eilmer
	- rm -rf geom
	- rm -rf $(ESTCN_DIR)/estcn-manual-for-hugo/estcn
	- rm -rf gas
	- rm -rf flow
	- rm -rf static/pdfs
	- rm -rf static/html
	- rm -rf public
	- rm $(EILMER_DIR)/eilmer-reference-manual-for-hugo.adoc
	- rm $(EILMER_DIR)/eilmer-catalog-for-hugo.adoc
	- rm $(EILMER_DIR)/eilmer-faq-for-hugo.adoc
	- rm $(EILMER_DIR)/technical-notes/hybrid-fluxes-for-hugo.adoc
	- rm -rf $(EILMER_DIR)/technical-notes/hybrid-fluxes-for-hugo
	- rm $(ESTCN_DIR)/estcn-manual-for-hugo.adoc
	- rm -rf $(ESTCN_DIR)/estcn-manual-for-hugo
	- rm $(L1D_DIR)/l1d-reference-manual-for-hugo.adoc
	- rm $(NENZF1D_DIR)/nenzf1d-manual-for-hugo.adoc
	- rm $(PYTHON_DIR)/gd-manual-for-hugo.adoc
	- rm $(PYTHON_DIR)/geometry-manual-for-hugo.adoc
	- rm $(PYTHON_DIR)/nm-manual-for-hugo.adoc


eilmer-catalog:
	cp ../doc/eilmer-catalog-of-examples.adoc $(EILMER_DIR)/
	mkdir -p $(EILMER_DIR)/eilmer-catalog-for-hugo
	cp -r ../doc/gallery $(EILMER_DIR)/eilmer-catalog-for-hugo/
	sed -i '/^:stylesheet:/d' $(EILMER_DIR)/eilmer-catalog-of-examples.adoc
	cat $(EILMER_DIR)/eilmer-catalog-hugo-header.txt $(EILMER_DIR)/eilmer-catalog-of-examples.adoc > $(EILMER_DIR)/eilmer-catalog-for-hugo.adoc
	rm $(EILMER_DIR)/eilmer-catalog-of-examples.adoc

eilmer-reference-manual:
	cp -r ../doc/eilmer/ .
	cp -r ../doc/geom/ .
	cp ../doc/eilmer-reference-manual.adoc $(EILMER_DIR)/
	sed -i '/^:stylesheet:/d' $(EILMER_DIR)/eilmer-reference-manual.adoc
	cat $(EILMER_DIR)/eilmer-reference-manual-hugo-header.txt $(EILMER_DIR)/eilmer-reference-manual.adoc > $(EILMER_DIR)/eilmer-reference-manual-for-hugo.adoc
	rm $(EILMER_DIR)/eilmer-reference-manual.adoc

eilmer-faq:
	cp ../doc/eilmer-faq.adoc $(EILMER_DIR)/
	sed -i '/^:stylesheet:/d' $(EILMER_DIR)/eilmer-faq.adoc
	cat $(EILMER_DIR)/eilmer-faq-hugo-header.txt $(EILMER_DIR)/eilmer-faq.adoc > $(EILMER_DIR)/eilmer-faq-for-hugo.adoc
	rm $(EILMER_DIR)/eilmer-faq.adoc

l1d-reference-manual:
	cp ../doc/l1d-reference-manual.adoc $(L1D_DIR)/
	sed -i '/^:stylesheet:/d' $(L1D_DIR)/l1d-reference-manual.adoc
	cat $(L1D_DIR)/l1d-reference-manual-hugo-header.txt $(L1D_DIR)/l1d-reference-manual.adoc > $(L1D_DIR)/l1d-reference-manual-for-hugo.adoc
	rm $(L1D_DIR)/l1d-reference-manual.adoc

estcn-manual:
	cp ../doc/estcn-manual.adoc $(ESTCN_DIR)/
	sed -i '/^:stylesheet:/d' $(ESTCN_DIR)/estcn-manual.adoc
	cp -r ../doc/estcn $(ESTCN_DIR)/estcn-manual-for-hugo/
	cat $(ESTCN_DIR)/estcn-manual-hugo-header.txt $(ESTCN_DIR)/estcn-manual.adoc > $(ESTCN_DIR)/estcn-manual-for-hugo.adoc
	rm $(ESTCN_DIR)/estcn-manual.adoc

nenzf1d-manual:
	cp ../doc/nenzf1d-manual.adoc $(NENZF1D_DIR)/
	sed -i '/^:stylesheet:/d' $(NENZF1D_DIR)/nenzf1d-manual.adoc
	cat $(NENZF1D_DIR)/nenzf1d-manual-hugo-header.txt $(NENZF1D_DIR)/nenzf1d-manual.adoc > $(NENZF1D_DIR)/nenzf1d-manual-for-hugo.adoc
	rm $(NENZF1D_DIR)/nenzf1d-manual.adoc

gas-dynamic-manual: geometry-manual nm-manual
	cp -r ../doc/gas/ .
	cp -r ../doc/flow/ .
	cp ../doc/library-reference-manual.adoc $(PYTHON_DIR)/
	sed -i '/^:stylesheet:/d' $(PYTHON_DIR)/library-reference-manual.adoc
	cat $(PYTHON_DIR)/gd-manual-hugo-header.txt $(PYTHON_DIR)/library-reference-manual.adoc > $(PYTHON_DIR)/gd-manual-for-hugo.adoc
	rm $(PYTHON_DIR)/library-reference-manual.adoc

geometry-manual:
	cp -r ../doc/geom/ .
	cp ../doc/geometry-reference-manual.adoc $(PYTHON_DIR)/
	sed -i '/^:stylesheet:/d' $(PYTHON_DIR)/geometry-reference-manual.adoc
	cat $(PYTHON_DIR)/geometry-manual-hugo-header.txt $(PYTHON_DIR)/geometry-reference-manual.adoc > $(PYTHON_DIR)/geometry-manual-for-hugo.adoc
	rm $(PYTHON_DIR)/geometry-reference-manual.adoc

nm-manual:
	cp ../doc/nm-reference-manual.adoc $(PYTHON_DIR)/
	sed -i '/^:stylesheet:/d' $(PYTHON_DIR)/nm-reference-manual.adoc
	cat $(PYTHON_DIR)/nm-manual-hugo-header.txt $(PYTHON_DIR)/nm-reference-manual.adoc > $(PYTHON_DIR)/nm-manual-for-hugo.adoc
	rm $(PYTHON_DIR)/nm-reference-manual.adoc

#--- Technical notes

hybrid-fluxes:
	cp -r ../doc/technical-notes/hybrid-fluxes.adoc $(EILMER_DIR)/technical-notes/
	- mkdir -p $(EILMER_DIR)/technical-notes/hybrid-fluxes-for-hugo
	cp -r ../doc/technical-notes/resources $(EILMER_DIR)/technical-notes/hybrid-fluxes-for-hugo
	sed -i '/^:stylesheet:/d' $(EILMER_DIR)/technical-notes/hybrid-fluxes.adoc
	cat $(EILMER_DIR)/technical-notes/hybrid-fluxes-hugo-header.txt $(EILMER_DIR)/technical-notes/hybrid-fluxes.adoc > $(EILMER_DIR)/technical-notes/hybrid-fluxes-for-hugo.adoc
	rm $(EILMER_DIR)/technical-notes/hybrid-fluxes.adoc



#--- User Guides

eilmer-user-guide:
	- mkdir -p static/pdfs
	cd $(DGD_DOC)/eilmer/user-guide; make pdf
	cp $(DGD_DOC)/eilmer/user-guide/eilmer-user-guide.pdf static/pdfs/

geometry-user-guide:
	- mkdir -p static/pdfs
	cd $(DGD_DOC)/geometry/user-guide; make pdf
	cp $(DGD_DOC)/geometry/user-guide/geometry-user-guide.pdf static/pdfs/

gas-user-guide:
	- mkdir -p static/pdfs
	cd $(DGD_DOC)/gas/user-guide; make pdf
	cp $(DGD_DOC)/gas/user-guide/gas-user-guide.pdf static/pdfs/

reacting-gas-guide:
	- mkdir -p static/pdfs
	cd $(DGD_DOC)/gas/reacting-gas; make pdf
	cp $(DGD_DOC)/gas/reacting-gas/reacting-gas-guide.pdf static/pdfs/

foammesh-user-guide:
	- mkdir -p static/pdfs
	cd $(DGD_DOC)/foammesh/user-guide; make pdf
	cp $(DGD_DOC)/foammesh/user-guide/foammesh-user-guide.pdf static/pdfs/

presentations:
	- mkdir -p static/pdfs
	rsync -av $(DGD_DOC)/presentations/2014-06-eilmer3-talk-pj/eilmer3-talk-june-2014.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2015-05-cfh-seminar-pj/eilmer4-talk-may-2015.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2015-11-accm-conference/eilmer4-talk-nov-2015.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2016-03-cfh-seminar-rjg/ static/html/
	rsync -av $(DGD_DOC)/presentations/2016-05-cfh-seminar-pj/eilmer4-talk-pj-may-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2016-09-paris-talk-pj/eilmer4-talk-pj-sep-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2016-12-oxford-talk-pj/eilmer4-talk-pj-dec-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2017-03-cfh-seminar-pj/eilmer4-talk-pj-march-2017.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2017-04-cfh-seminar-rjg/seminar-2017-04-06.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2017-04-cfh-seminar-pj/eilmer4-talk-pj-april-2017.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-05-nus-short-course-rjg/tex/nus-short-course-on-eilmer-day-1.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-05-nus-short-course-rjg/tex/nus-short-course-on-eilmer-day-2.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-06-strathclyde-talk-rjg/tex/rjg-seminar-strathclyde-2018-06-08.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-07-cfh-seminar-pj/eilmer-talk-pj-july-2018.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-07-cfh-seminar-rjg/tex/cfh-seminar-jul-2018.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2019-08-cfh-seminar-pj/eilmer-talk-pj-aug-2019.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2019-10-03-cfh-seminar-rjg/tex/cfh-seminar-oct-2019.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2011-shock-tube-workshop/nenzfr-presentation-pj-iwstt-2011.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2020-12-cfh-seminar-pj/eilmer-talk-pj-2020-dec.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2021-01-cfh-seminar-rjg/cfh-seminar-jan-2021.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2021-05-cfh-seminar-kd/presentation.pdf static/pdfs/kad-cfh-seminar-05-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-05-cfh-seminar-nng/talk.pdf static/pdfs/nng-cfh-seminar-05-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-07-cfh-seminar-cmj/james-cfh-talk-2021.pdf static/pdfs/cmj-cfh-talk-07-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-07-uqtum-seminar-nng/talk.pdf static/pdfs/nng-uqtum-seminar-07-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-07-uqtum-seminar-kd/talk.pdf static/pdfs/kd-uqtum-seminar-07-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-09-cfh-seminar-viv-bone-andrew-lock/main.pdf static/pdfs/viv-bone-andrew-lock-cfh-seminar-09-2021.pdf
	rsync -av $(DGD_DOC)/presentations/2021-11-eilmer-meetup-fabs/eilmer_moving_grids.pdf static/pdfs/fabs-flying-cubes-eilmer-meetup-2021-11-25.pdf
	rsync -av $(DGD_DOC)/presentations/2021-11-3rd-franco-oz-symposium-rjg/third-franco-oz-hyp-symp-2021-11-25-rjg.pdf static/pdfs/third-franco-oz-hyp-symp-2021-11-25-rjg.pdf

theses:
	- mkdir -p static/pdfs
	rsync -av $(DGD_DOC)/theses/adriaan-window-MPhil-thesis-oct-2008.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/andrew-denman-PhD-thesis-jan-2007.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/ben-stewart-phd-thesis-oct-2004.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/brendan-oflaherty-phd-thesis-june-2012.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/carolyn-jacobs-phd-finalthesis-UQversion-aug-2011.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/charles-lilley-phd-thesis-jun-2005.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/chris-craddock-phd-thesis-aug-1999.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/dan-potter-phd-thesis-may-2011.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/david-gildfind-phd-thesis-oct-2012.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/dwishen-ramanah-mphil_thesis-jan-2007.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/elise-fahy-phd-thesis-dec-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/haruko-ishikawa-phd-thesis-dec-1999.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/ian-johnston-honours-thesis-1994.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/ian-johnston-phd-thesis-jan-1999-export.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/james-faddy-masters-thesis-aug-2000.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/jens-kunze-phd-thesis-feb-2020.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/joseph-tang-PhD-thesis-jun-2008.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/kan-qin-phd-thesis-oct-2017.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/kevin-austin-phd-thesis-june-2002.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/kyle-damm-phd-thesis-feb-2020.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/matt-mcgilvray-PhD-thesis-sep-2008.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/matthew-smith-PhD-thesis-june-2008.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/michael-elford-masters-thesis-sep-2005.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/michael-scott-phd-thesis-june-2006.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/n33694032_PhD_submission.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/neil-mudford-phd-thesis-production-of-pulsed-nozzle-flows.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/paul-petrie-repar-phd-thesis-dec-1997.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/peter-jacobs-phd-thesis-may-1987.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/pierpaolo-toniato-phd-thesis-2019.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/richard-goozee-phd-thesis-apr-2003.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/rowan-gollan-PhD-thesis-feb-2009.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/vince-wheatley-masters-thesis-aug-2001.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/wilson-chan-phd-thesis-aug-2012.pdf static/pdfs/
