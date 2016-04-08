all: pooled_error.pdf suppmat.pdf

pooled_error.pdf: pooled_error.tex main_failsim.pdf main_okaysim.pdf main_power.pdf quantile_with.pdf main_power_zoom.pdf references.bib
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

main_power.pdf: ../simulations/ESpresso/results_power/ROC_1.eps
	epstopdf $< -o=$@

main_power_zoom.pdf: ../simulations/ESpresso/results_power/ROC_zoom_1.eps
	epstopdf $< -o=$@

extra_power_2.pdf: ../simulations/ESpresso/results_power/ROC_2.eps
	epstopdf $< -o=$@

extra_power_3.pdf: ../simulations/ESpresso/results_power/ROC_3.eps
	epstopdf $< -o=$@

extra_power_4.pdf: ../simulations/ESpresso/results_power/ROC_4.eps
	epstopdf $< -o=$@

extra_FDR.pdf: ../simulations/ESpresso/results_power/FDR.eps
	epstopdf $< -o=$@

quantile_with.pdf: ../quantile/with.eps
	epstopdf $< -o=$@

suppmat.pdf: suppmat.tex pooled_error.pdf extra_failsim.pdf extra_okaysim.pdf extra_power_2.pdf extra_power_3.pdf extra_power_4.pdf extra_FDR.pdf \
			quantile_schematic.pdf realdata_top_edgeR_ola1_sum.pdf realdata_top_edgeR_ola2_sum.pdf
	pdflatex suppmat.tex

quantile_schematic.pdf: ../quantile/sim_schematic.pdf
	cp $< $@

realdata_top_edgeR_ola1_sum.pdf: ../realdata/top_edgeR_ola1_sum.pdf
	cp $< $@

realdata_top_edgeR_ola2_sum.pdf: ../realdata/top_edgeR_ola2_sum.pdf
	cp $< $@
