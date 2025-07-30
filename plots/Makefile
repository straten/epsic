
FIGS := polarization_ellipse-crop.pdf polarization_sphere-crop.pdf polarization_cases-crop.pdf

all: $(FIGS)

%-crop.pdf: %.pdf
	pdfcrop $< $@

%.pdf: %.py
	python3 $< $@ 

