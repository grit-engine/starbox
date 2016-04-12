# (c) David Cunningham 2011, Licensed under the MIT license: http://www.opensource.org/licenses/mit-license.php

starbox: starbox.c
	gcc -std=c99 -pedantic -O2 -finline-functions -fomit-frame-pointer -march=prescott -Wall -Wextra $(CFLAGS) $< -lm $(LDFLAGS) -o $@

pnms: starbox starfile.csv
	./starbox -f 1024 -d 23.5 -r 90 -g 1 -S 20 starfile.csv
	touch pnms

%.pnm: %.png
	pngtopnm $< > $@

%.pgm: %.png
	pngtopnm -alpha $< > $@

%u.pnm: %.pnm
	pnmflip -r180 $< > $@

TMP1:=/tmp/wtflol1.txt
TMP2:=/tmp/wtflol2.txt

starfield.png: pnms starfield_template.pgm starfield_template.pnm
	$(MAKE) face5u.pnm
	#$(MAKE) face0u.pnm face1u.pnm face2u.pnm face3u.pnm face4u.pnm face5.pnm
	cat starfield_template.pnm > $(TMP1)
	cat $(TMP1) | pnmcomp -xoff 3072 -yoff    0 face0.pnm > $(TMP2)
	cat $(TMP2) | pnmcomp -xoff 2048 -yoff    0 face1.pnm > $(TMP1)
	cat $(TMP1) | pnmcomp -xoff 1024 -yoff    0 face2.pnm > $(TMP2)
	cat $(TMP2) | pnmcomp -xoff    0 -yoff    0 face3.pnm > $(TMP1)
	cat $(TMP1) | pnmcomp -xoff 3072 -yoff 1024 face4.pnm > $(TMP2)
	cat $(TMP2) | pnmcomp -xoff 1024 -yoff 1024 face5u.pnm > $(TMP1)
	pnmtopng -alpha starfield_template.pgm $(TMP1) > $@
	rm $(TMP1) $(TMP2)



clean:
	rm -f *.pnm *.pgm starfield.png test

squeakyclean: clean
	rm -f starbox

.PHONY: clean squeakyclean

# vim: ts=8:sw=8:noet
