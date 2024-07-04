prep-gas:
	prep-gas ideal-air.inp ideal-air-gas-model.lua

prep-puffin:
	puffin-prep --job=busemann-diffuser

prep: prep-gas prep-puffin

run:
	puffin --job=busemann-diffuser

vtk:
	puffin-post --job=busemann-diffuser --output=vtk

wall-data:
	puffin-post --job=busemann-diffuser --output=stream --cell-index=$ --stream-index=0
