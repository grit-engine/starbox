starbox: starbox.c
	gcc -std=c99 -pedantic -O2 -finline-functions -fomit-frame-pointer -march=prescott $(CFLAGS) $< -lm $(LDFLAGS) -o $@

test:
	./starbox starfile.csv -f 1024 -d 23.5 -r -90 -g 3 -S 8

%.png: %.pnm
	pnmtopng $< > $@

pngs:
	$(MAKE) test
	$(MAKE) face0.png face1.png face2.png face3.png face4.png face5.png

clean:
	rm -f *.png *.pnm

squeakyclean: clean
	rm -f starbox

.PHONY: test pngs clean squeakyclean

# vim: ts=8:sw=8:noet
