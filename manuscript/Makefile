all: pooled_error.pdf suppmat.pdf

pooled_error.pdf: pooled_error.tex main_failsim.pdf main_okaysim.pdf main_power.pdf main_power_zoom.pdf main_swapped.pdf references.bib
	pdflatex pooled_error.tex
	bibtex pooled_error
	pdflatex pooled_error.tex
	pdflatex pooled_error.tex 

main_failsim.pdf: ../simulations/ESpresso/results_failsim/main_raw.eps
	epstopdf $< -o=$@

main_okaysim.pdf: ../simulations/ESpresso/results_failsim/main_sum.eps
	epstopdf $< -o=$@

main_power.pdf: ../simulations/ESpresso/results_power/ROC_1.eps
	epstopdf $< -o=$@

main_power_zoom.pdf: ../simulations/ESpresso/results_power/ROC_zoom_1.eps
	epstopdf $< -o=$@

main_swapped.pdf: ../realdata/all_scrambled.eps
	epstopdf $< -o=$@

suppmat.pdf: suppmat.tex pooled_error.pdf \
		../simulations/ESpresso/results_failsim/supp_raw.pdf ../simulations/ESpresso/results_failsim/supp_sum.pdf \
		../simulations/ESpresso/results_power/ROC_2.pdf ../simulations/ESpresso/results_power/ROC_3.pdf ../simulations/ESpresso/results_power/ROC_4.pdf ../simulations/ESpresso/results_power/FDR.pdf \
		../realdata/ESpresso/real_ranks.pdf ../realdata/ESpresso/top_edgeR_3_sum.pdf ../realdata/ESpresso/top_edgeR_4_sum.pdf
	pdflatex suppmat.tex
	pdflatex suppmat.tex

