
help:
	@echo "Main targets are:"
	@echo -e "   demo"
	@echo -e "   clean"

run-pc:
	lmr-verify --norms="rho" --case-tag="predictor_corrector"


run-rk3:
	lmr-verify --norms="rho" --case-tag="classic_rk3"

run-rk4:
	lmr-verify --norms="rho" --case-tag="classic_rk4"

run-rkl2:
	lmr-verify --norms="rho" --case-tag="rkl2"

ooa-pc:
	@echo -n "Observed Order of Accuracy (Predictor Corrector): "
	@cat case-predictor_corrector/observed-order-case-predictor_corrector.dat | awk 'NR == 2 {print $$3}'

ooa-rk3:
	@echo -n "Observed Order of Accuracy (Classic RK3):         "
	@cat case-classic_rk3/observed-order-case-classic_rk3.dat | awk 'NR == 2 {print $$3}'

ooa-rk4:
	@echo -n "Observed Order of Accuracy (Classic RK4):         "
	@cat case-classic_rk4/observed-order-case-classic_rk4.dat | awk 'NR == 2 {print $$3}'

ooa-rkl2:
	@echo -n "Observed Order of Accuracy (STS RKL2):            "
	@cat case-rkl2/observed-order-case-rkl2.dat | awk 'NR == 2 {print $$3}'

demo: run-pc run-rk3 run-rk4 run-rkl2 ooa-pc ooa-rk3 ooa-rk4 ooa-rkl2
	@echo "Done"

clean:
	- rm -rf case-*

