cheatsheet: eilmer-cheatsheet-lua.pdf eilmer-cheatsheet-cmds.pdf
	qpdf --empty --pages eilmer-cheatsheet-cmds.pdf eilmer-cheatsheet-lua.pdf -- eilmer-cheatsheet-v5.0.0.pdf

eilmer-cheatsheet-lua.pdf:
	typst c --font-path=../../../gdtk-docs/fonts eilmer-cheatsheet-lua.typ

eilmer-cheatsheet-cmds.pdf:
	typst c --font-path=../../../gdtk-docs/fonts eilmer-cheatsheet-cmds.typ

clean:
	- rm eilmer-cheatsheet-*.pdf
