all: pooled_error.pdf suppmat.pdf

pooled_error.pdf: pooled_error.tex main_failsim.pdf main_okaysim.pdf main_power.pdf main_power_zoom.pdf main_swapped.pdf references.bib
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

main_swapped.pdf: ../realdata/all_scrambled.eps
	epstopdf $< -o=$@

extra_rankings.pdf: ../realdata/ESpresso/real_ranks.eps
	epstopdf $< -o=$@

suppmat.pdf: suppmat.tex pooled_error.pdf extra_failsim.pdf extra_okaysim.pdf extra_power_2.pdf extra_power_3.pdf extra_power_4.pdf extra_FDR.pdf \
		extra_rankings.pdf realdata_top_edgeR_ola1_sum.pdf realdata_top_edgeR_ola2_sum.pdf
	pdflatex suppmat.tex

quantile_schematic.pdf: ../quantile/sim_schematic.pdf
	cp $< $@

realdata_top_edgeR_ola1_sum.pdf: ../realdata/top_edgeR_ola1_sum.pdf
	cp $< $@

realdata_top_edgeR_ola2_sum.pdf: ../realdata/top_edgeR_ola2_sum.pdf
	cp $< $@
