all: pooled_error.pdf suppmat.pdf

pooled_error.pdf: pooled_error.tex main_failsim.pdf main_okaysim.pdf power_with.pdf quantile_with.pdf power_with_zoom.pdf references.bib
	pdflatex pooled_error.tex
	bibtex pooled_error
	pdflatex pooled_error.tex
	pdflatex pooled_error.tex 

main_failsim.pdf: ../simulations/ESpresso/results_failsim/main_raw.eps
	epstopdf $< -o=$@

extra_failsim.pdf: ../simulations/ESpresso/results_failsim/supp_raw.eps
	epstopdf $< -o=$@

main_okaysim.pdf: ../simulations/ESpresso/results_failsim/main_sum.eps
	epstopdf $< -o=$@

extra_okaysim.pdf: ../simulations/ESpresso/results_failsim/supp_sum.eps
	epstopdf $< -o=$@

power_with.pdf: ../power/with.eps
	epstopdf $< -o=$@

power_with_zoom.pdf: ../power/with_zoom.eps
	epstopdf $< -o=$@

quantile_with.pdf: ../quantile/with.eps
	epstopdf $< -o=$@

suppmat.pdf: suppmat.tex pooled_error.pdf extra_failsim.pdf extra_okaysim.pdf quantile_schematic.pdf realdata_top_edgeR_ola1_sum.pdf realdata_top_edgeR_ola2_sum.pdf
	pdflatex suppmat.tex

quantile_schematic.pdf: ../quantile/sim_schematic.pdf
	cp $< $@

realdata_top_edgeR_ola1_sum.pdf: ../realdata/top_edgeR_ola1_sum.pdf
	cp $< $@

realdata_top_edgeR_ola2_sum.pdf: ../realdata/top_edgeR_ola2_sum.pdf
	cp $< $@
