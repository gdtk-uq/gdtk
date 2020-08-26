# makefile for the new web site, 2020-05-05.

DGD_DOC ?= ../../dgd-docs
PDFS := eilmer-user-guide \
	geometry-user-guide \
	gas-user-guide \
	reacting-gas-guide \
	foammesh-user-guide \
	presentations \
	theses

HTML := eilmer-reference-manual.html \
	l1d-reference-manual.html \
	library-reference-manual.html \
	geometry-reference-manual.html \
	nm-reference-manual.html \
	estcn-manual.html \
	build-lut-manual.html \
	eilmer-catalog-of-examples.html \
	eilmer-faq.html

build: $(PDFS) $(HTML) images
	hugo

sync:
	rsync -zz -vrlt -p --chmod=Dg+s,ug+w,Fo-w,+X -g --groupmap=*:other public/ cfcfd:/opt/local/share/httpd/htdocs/

clean:
	cd ../doc; make clean
	- rm -r static
	- rm -r public

images:
	- mkdir static/images
	cp ../images/* static/images/
	cp ../images/eilmer-icon-2.png static/favicon.png

eilmer-reference-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

eilmer-catalog-of-examples.html:
	- mkdir -p static/html
	- mkdir -p static/html/gallery
	cd ../doc; make $@
	cp ../doc/$@ static/html/
	cp ../doc/gallery/*.png static/html/gallery/

eilmer-faq.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

l1d-reference-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

estcn-manual.html:
	- mkdir -p static/html
	- mkdir -p static/html/estcn
	cd ../doc; make $@
	cp ../doc/$@ static/html/
	cp ../doc/estcn/*.png static/html/estcn/

build-lut-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

library-reference-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

geometry-reference-manual.html:
	- mkdir -p static/html
	- mkdir -p static/html/geom
	cd ../doc; make $@
	cp ../doc/$@ static/html/
	cp ../doc/geom/*.png static/html/geom/

nm-reference-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

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
	rsync -av $(DGD_DOC)/presentations/2016-05-cfh-seminar-pj/eilmer4-talk-pj-may-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2016-09-paris-talk-pj/eilmer4-talk-pj-sep-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2016-12-oxford-talk-pj/eilmer4-talk-pj-dec-2016.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2017-03-cfh-seminar-pj/eilmer4-talk-pj-march-2017.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2017-04-cfh-seminar-pj/eilmer4-talk-pj-april-2017.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-05-nus-short-course-rjg/tex/nus-short-course-on-eilmer-day-1.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-05-nus-short-course-rjg/tex/nus-short-course-on-eilmer-day-2.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-06-strathclyde-talk-rjg/tex/rjg-seminar-strathclyde-2018-06-08.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-07-cfh-seminar-pj/eilmer-talk-pj-july-2018.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2018-07-cfh-seminar-rjg/tex/cfh-seminar-jul-2018.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2019-08-cfh-seminar-pj/eilmer-talk-pj-aug-2019.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2019-10-03-cfh-seminar-rjg/tex/cfh-seminar-oct-2019.pdf static/pdfs/
	rsync -av $(DGD_DOC)/presentations/2011-shock-tube-workshop/nenzfr-presentation-pj-iwstt-2011.pdf static/pdfs/

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
	rsync -av $(DGD_DOC)/theses/joseph-tang-PhD-thesis-jun-2008.pdf static/pdfs/
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
	rsync -av $(DGD_DOC)/theses/richard-goozee-phd-thesis-apr-2003.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/rowan-gollan-PhD-thesis-feb-2009.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/vince-wheatley-masters-thesis-aug-2001.pdf static/pdfs/
	rsync -av $(DGD_DOC)/theses/wilson-chan-phd-thesis-aug-2012.pdf static/pdfs/
