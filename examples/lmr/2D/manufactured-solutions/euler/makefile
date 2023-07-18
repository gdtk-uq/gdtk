
help:
	@echo "Main targets are:"
	@echo -e "   demo"
	@echo -e "   plot"
	@echo -e "   clean"

demo:
	python verify.py --norms="rho,T,p,vel.x,vel.y" --number-cells="8,16,32"

plot:
	python verify.py --case-number=0 --norms="rho,T,p,vel.x,vel.y" --number-cells="8,16,32,64,128,256"
	gnuplot -c plot-norms-n-obs-order.gplot case-000

clean:
	- rm -rf case-000

