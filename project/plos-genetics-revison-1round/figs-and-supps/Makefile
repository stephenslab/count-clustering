
supps: supplementals.tex ../header.tex S1_Fig.tex S2_Fig.tex S3_Fig.tex S4_Fig.tex S5_Fig.tex S1_Table.tex S2_Table.tex S3_Table.tex
	pdflatex supplementals.tex

fig1: fig1.tex ../header.tex ../../plots/gtex-figures/gtex-k20-05-01-2016.pdf
	pdflatex fig1.tex

fig2: fig2.tex ../header.tex ../../plots/hierarchy_F_thin_0_1.png ../../plots/admixture_F_thin_0_1.png
	pdflatex fig2.tex

fig3: fig3.tex ../header.tex ../../plots/jaitin-figures/jaitin-main2.pdf
	pdflatex fig3.tex

fig4: fig4.tex ../header.tex ../../plots/deng-figures/deng-main-v4.pdf
	pdflatex fig4.tex

figS1: S1_Fig.tex ../header.tex ../../plots/gtex-figures/gtex-multiple-ks-04-30-2016.pdf
		pdflatex S1_Fig.tex

figS2: S2_Fig.tex ../header.tex ../../plots/gtex-figures/gtex_thinned_04_25_2016.pdf
		pdflatex S2_Fig.tex

figS3: S3_Fig.tex ../header.tex ../../plots/rsz_1hierarchy_F_thin_0_01.png ../../plots/rsz_1hierarchy_F_thin_0_01.png ../../plots/rsz_1hierarchy_F_thin_0_0001.png ../../plots/rsz_1admixture_F_thin_0_0001.png
		pdflatex S3_Fig.tex

figS4: S4_Fig.tex ../header.tex ../../src/figure/gtex-pca-vs-tsne.Rmd/gtex-with-legend.jpeg
		pdflatex S4_Fig.tex

figS5: S5_Fig.tex ../header.tex ../../src/figure/deng-pca-vs-tsne.Rmd/deng-with-legend.jpeg
		pdflatex S5_Fig.tex

tableS1: S1_Table.tex ../header.tex
		pdflatex S1_Table.tex

tableS2: S2_Table.tex ../header.tex
		pdflatex S2_Table.tex

tableS3: S3_Table.tex ../header.tex
		pdflatex S3_Table.tex

clean:
	\rm -f *.aux *.bbl *.blg *.log *.bak *.out *~ *.Rout */*~ */*.Rout */*.aux */*.log

cleanall:
	\rm -f *.aux *.bbl *.blg *.log *.bak *.out *~ *.Rout *.pdf */*~ */*.Rout */*.aux */*.log
