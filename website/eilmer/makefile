DGD_DOC ?= ../../../dgd-doc
PDFS := eilmer-user-guide \
	geometry-user-guide \
	gas-user-guide \
	foammesh-user-guide

HTML := users-reference-manual.html

build: $(PDFS) $(HTML) images
	hugo

sync:
	rsync -zav public/ cfcfd:/opt/local/share/httpd/htdocs/eilmer/

images:
	cp ../../images/* static/images/

users-reference-manual.html:
	cd ../../doc; make $@
	cp ../../doc/$@ static/html/

eilmer-user-guide:
	cd $(DGD_DOC)/eilmer/user-guide; make pdf
	cp $(DGD_DOC)/eilmer/user-guide/eilmer-user-guide.pdf static/pdfs/

geometry-user-guide:
	cd $(DGD_DOC)/geometry/user-guide; make pdf
	cp $(DGD_DOC)/geometry/user-guide/geometry-user-guide.pdf static/pdfs/

gas-user-guide:
	cd $(DGD_DOC)/gas/user-guide; make pdf
	cp $(DGD_DOC)/gas/user-guide/gas-user-guide.pdf static/pdfs/

foammesh-user-guide:
	cd $(DGD_DOC)/foammesh/user-guide; make pdf
	cp $(DGD_DOC)/foammesh/user-guide/foammesh-user-guide.pdf static/pdfs/
