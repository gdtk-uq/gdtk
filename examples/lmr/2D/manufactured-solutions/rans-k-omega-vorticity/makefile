help:
	@echo "Main targets are:"
	@echo -e "   demo"
	@echo -e "   plot"
	@echo -e "   clean"

demo:
	lmr-verify --norms="rho,T,p,vel.x,vel.y,tq-tke,tq-omega" --grid-levels="5,4"

verify:
	lmr-verify --case-number=0 --norms="rho,T,p,vel.x,vel.y,tq-tke,tq-omega" 

plot:
	gnuplot -c plot-norms-n-obs-order.gplot case-000

clean:
	- rm -rf case-000

