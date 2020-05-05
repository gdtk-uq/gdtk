# makefile for the new web site, 2020-05-05.

DGD_DOC ?= ../../dgd-docs
PDFS := eilmer-user-guide \
	geometry-user-guide \
	gas-user-guide \
	reacting-gas-guide \
	foammesh-user-guide \
	presentations

HTML := users-reference-manual.html l1d-reference-manual.html

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

users-reference-manual.html:
	- mkdir -p static/html
	cd ../doc; make $@
	cp ../doc/$@ static/html/

l1d-reference-manual.html:
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
