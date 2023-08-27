help:
	@echo "Main targets are:"
	@echo -e "   demo"
	@echo -e "   plot"
	@echo -e "   clean"

demo:
	lmr-verify --norms="rho,T,p,vel.x,vel.y,vel.z" --grid-levels="4,3"

verify:
	lmr-verify --case-number=0 --norms="rho,T,p,vel.x,vel.y,vel.z" 

plot:
	gnuplot -c plot-norms-n-obs-order.gplot case-000

clean:
	- rm -rf case-000

