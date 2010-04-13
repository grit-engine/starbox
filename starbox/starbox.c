/* Copyright (c) Graeme Cole 2010
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */



#define _GNU_SOURCE

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#define DEFAULT_FACE_DIM 400
#define DEFAULT_MID_MAG 3
#define DEFAULT_SCALE_FACTOR 10
#define DEFAULT_CUTOFF_MAG 10000

#define OFMT_PNM 1
#define OFMT_C 2

#define FACE_FRONT 0
#define FACE_LEFT 1
#define FACE_BACK 2
#define FACE_RIGHT 3
#define FACE_TOP 4
#define FACE_BOTTOM 5

#define MAG_LOG_BASE 2.512


int verbose = 0;


double pixel_brightness(double mag, double mid) {
    double brightness;

    /* An N-magnitude star is MAG_LOG_BASE times brighter than an
       (N+1)-magnitude star. */

    if (mag <= mid) {
        double multiplier = pow(MAG_LOG_BASE, mid - mag);
        brightness = 0.5 * multiplier;
    }
    else {
        double divisor = pow(MAG_LOG_BASE, mag - mid);
        brightness = 0.5 / divisor;
    }

    if (brightness < 0)
        brightness = 0;
    if (brightness > 1)
        brightness = 1;

    return brightness;
}

void ra_dec_to_xyz_aux(double ra, double dec, double *x, double *y, double *z, int magnify) {
    double cosdec;
    double normfac;

    *z = sin(dec * M_PI / 180);

    cosdec = cos(dec * M_PI / 180);
    *x = cosdec * cos(ra * M_PI / 180);
    *y = cosdec * sin(ra * M_PI / 180);

    if (magnify) {
        /* Normalise these vectors so that the one with the highest
           magnitude is a unit vector */
        
        if (fabs(*x) > fabs(*y)) {
            if (fabs(*x) > fabs(*z)) {
                normfac = 1 / fabs(*x);
            }
            else {
                normfac = 1 / fabs(*z);
            }
        }
        else if (fabs(*y) > fabs(*z)) {
            normfac = 1 / fabs(*y);
        }
        else {
            normfac = 1 / fabs(*z);
        }

        *x *= normfac;
        *y *= normfac;
        *z *= normfac;
    }
}

void ra_dec_to_xyz(double ra, double dec, double *x, double *y, double *z) {
    ra_dec_to_xyz_aux(ra, dec, x, y, z, 0);
}

void ra_dec_to_xyz_magnify(double ra, double dec, double *x, double *y, double *z) {
    ra_dec_to_xyz_aux(ra, dec, x, y, z, 1);
}

void xyz_to_ra_dec(double x, double y, double z, double *ra, double *dec) {
    double w;

    if (x == 0) {
        if (y < 0)
            *ra = 3 * M_PI / 2;
        else
            *ra = M_PI / 2;
    }
    else if (y == 0) {
        *ra = 0;
    }
    else {
        *ra = atan(y / x) * 180.0 / M_PI;
        if (x < 0) {
            *ra += 180;
        }
        if (*ra >= 360)
            *ra -= 360;
        if (*ra < 0)
            *ra += 360;
    }

    if (z == 0) {
        *dec = 0;
    }
    else {
        w = sqrt(x * x + y * y);

        if (w == 0) {
            *dec = (z > 0) ? 90 : -90;
        }
        else {
            *dec = atan(fabs(z) / w) * 180.0 / M_PI;
            if (z < 0)
                *dec = -*dec;
        }
    }
}


/* Read a line from the CSV file. The star's right-ascension should be in the
   range [0,24), the declination [-90, 90].
   Return 1 on success, -1 for a malformed record or 0 for end of file.
*/
int
read_star_line(FILE *file, double *ra, double *dec, double *mag, int line_no,
        int ra_field, int dec_field, int mag_field) {
    double ra_24;
    int ret = 0;
    int c;

    int field = 0;

    c = 0;
    while (c != '\n' && c != EOF) {
        if (c == ',' || c == 0) {
            ++field;
            /* Right Ascension, between 0 and 24 */
            if (field == ra_field) {
                ret = fscanf(file, " %lf", &ra_24);
            }
            /* Declination, between -90 and 90 */
            else if (field == dec_field) {
                ret = fscanf(file, " %lf", dec);
            }
            /* Magnitude */
            else if (field == mag_field) {
                ret = fscanf(file, " %lf", mag);
            }

            /* If we read the field incorrectly, it's a malformed
               line. If ret < 0, then we've reached EOF. */
            if (ret == 0) {
                goto malformed;
            }
            else if (ret < 0) {
                c = EOF;
                continue;
            }
        }
        c = fgetc(file);
    }

    if (c == EOF)
        return 0;

    /* Convert right ascension from 0-24 to 0-360 */
    *ra = ra_24 * 15;

    if (*ra < 0 || *ra >= 360) {
        error(0, 0, "line %d: invalid right ascension %f (must be 0 <= x < 24)\n", line_no, ra_24);
        return -1;
    }
    if (*dec < -90 || *dec > 90) {
        error(0, 0, "line %d: invalid declination %f (must be -90 <= x <= 90)\n", line_no, *dec);
        return -1;
    }

    return 1;

malformed:
    error(0, 0, "line %d malformed.\n", line_no);

    /* Read to the next newline */
    while ((c = fgetc(file)) != '\n' && c != EOF);

    return -1;
}

void
progress_report(int stars_processed, int stars_drawn, int stars_discarded,
        int malformed_lines, int removed_sun, int clipped,
        int final_report) {
    fprintf(stderr, "%7d star%s read  %7d drawn",
            stars_processed, stars_processed == 1 ? "" : "s",
            stars_drawn);
    if (stars_discarded > 0)
        fprintf(stderr, "  %7d discarded", stars_discarded);
    if (malformed_lines > 0)
        fprintf(stderr, "  %7d record%s rejected", malformed_lines,
                malformed_lines == 1 ? "" : "s");
    if (final_report) {
        if (removed_sun)
            fprintf(stderr, "     removed sun");
        fprintf(stderr, "\n");
        if (clipped > 0)
            fprintf(stderr, "%7d pixel%s clipped to maximum brightness\n", clipped, clipped == 1 ? "" : "s");
    }
    else
        fprintf(stderr, "\r");
}


void draw_star(float **faces, double ra, double dec, double mag,
        double ra_delta, double dec_delta, int face_dim,
        double mid_mag, int *clipped) {
    int face;
    float *pixel_dest;
    double x, y, z;
    double hvec, vvec;
    int image_x, image_y;

    /* Apply right ascension adjustment */
    ra += ra_delta;
    while (ra >= 360)
        ra -= 360;
    while (ra < 0)
        ra += 360;

    /* Apply declination adjustment, which is a bit more tricky than just
       adding a value.
       If the observer, facing RA 0, adjusts their declination, e.g. leans
       back 30 degrees, then stars at RA 0 just have their declination
       adjusted accordingly. Stars at 90 and 270 degrees RA stay exactly
       where they are, but other stars move dec_delta degrees around a
       circle whose centre axis runs between RA 270 and RA 90, that is,
       through the left and right sides of the observer. Basically this
       means the star's apparent right ascension will change as well as its
       declination. */

    if (dec_delta != 0) {
        if (ra == 0) {
            dec += dec_delta;
        }
        else if (ra == 180) {
            dec -= dec_delta;
        }
        else {
            /* Convert ra,dec to x,y,z */
            double x, y, z;
            double circ_angle;
            double hdist, radius;

            ra_dec_to_xyz(ra, dec, &x, &y, &z);

            /* If we call the radius of our sphere
               of stars "1", then the radius of the circle around
               which this star moves is determined by the star's
               declination and right ascension; a star atop the
               observer or directly in front would move round a
               larger circle than one that's to the right of the
               observer and near the horizon.
               
               The angle between the X vector and the new position
               of this star on the circle will be circ_angle.
             */


            /*
               The shortest horizontal distance between the
               observer and the centre of this circle is given by
               the Y value if you convert (ra, dec) into (x,y,z)
               where we'll call the distance to the star "1".  This
               is given by cos(dec)sin(ra).  Call this hdist.

               The angle between the horizontal through the
               observer and the edge of the circle is
               arccos(hdist).

               So the radius of our circle is sin(arccos(hdist)).
             */


            hdist = cos(dec * M_PI / 180.0) * sin(ra * M_PI / 180.0);

            assert(hdist >= -1);
            assert(hdist <= 1);

            radius = sin(acos(hdist));

            if (verbose)
                printf("hdist=%f, radius=%f\n", hdist, radius);

            if (radius != 0) {
                if (fabs(z) > fabs(radius)) {
                    /* Have seen fabs(z) > fabs(radius)
                       by a teeny tiny amount, probably
                       due to floating point pain */
                    z = fabs(radius) * (z > 0 ? 1 : -1);
                }
                assert(fabs(z) <= fabs(radius));
                circ_angle = asin(fabs(z) / fabs(radius)) * 180 / M_PI;
                if (x < 0) {
                    if (circ_angle > 0)
                        circ_angle = 180 - circ_angle;
                    else
                        circ_angle = -180 - circ_angle;
                }
                if (z < 0)
                    circ_angle = -circ_angle;
                circ_angle += dec_delta;

                if (circ_angle >= 360)
                    circ_angle -= 360;
                if (circ_angle < 0)
                    circ_angle += 360;

                /* Modify the star's X and Z co-ordinates.
                   The Y co-ord stays the same. */
                x = radius * cos(circ_angle * M_PI / 180.0);
                z = radius * sin(circ_angle * M_PI / 180.0);

                /* Convert back to RA/dec */
                xyz_to_ra_dec(x, y, z, &ra, &dec);
                
                if (verbose)
                    printf("x=%f, y=%f, z=%f\n", x, y, z);
            }
        }

        /* Make sure RA is between 0 and 360 */
        while (ra >= 360)
            ra -= 360;
        while (ra < 0)
            ra += 360;
        
        /* Make sure declination is between -90 and 90, and if we
           (e.g.) change -100 into -80 remember to flip the right
           ascension direction as well. */

        if (dec > 90) {
            dec = 180 - dec;
            ra = ra + 180;
            if (ra >= 360)
                ra -= 360;
        }
        if (dec < -90) {
            dec = -180 - dec;
            ra = ra + 180;
            if (ra >= 360)
                ra -= 360;
        }
    }


    if (verbose)
        printf("ra=%f,dec=%f,mag=%f\n", ra, dec, mag);

    /* Convert RA and Dec to normalised x,y,z, where the
       largest component is -1 or +1. */
    ra_dec_to_xyz_magnify(ra, dec, &x, &y, &z);

    /* Which face does this hit? */

    /* Face 0's centre is RA 0, Dec 0 (vector (1, 0, 0))
       Face 1's centre is RA 90, Dec 0 (vector (0, 1, 0)).
       Face 2's centre is RA 180, Dec 0 (vector (-1, 0, 0)).
       Face 3's centre is RA 270, Dec 0 (vector (0, -1, 0)).
       Face 4's centre is Dec 90 (directly above observer).
       Top row of the image is behind the observer if
       observer faces towards RA 0 (vector (0, 0, 1)).
       Face 5's centre is Dec -90 (directly below observer).
       Top row of the image is in front of the observer if
       observer faces towards RA 0 (vector (0, 0, -1)).
     */

    /* Work out which face we're drawing on. This is
       determined by which of x, y and z has the largest
       magnitude. */

    /* Set hvec to the vector component which will affect
       which column on the face we light up, and vvec to
       the vector component which will affect which row
       on the face we light up. */
    if (fabs(x) > fabs(y)) {
        if (fabs(x) > fabs(z)) {
            if (x > 0) {
                face = FACE_FRONT;
                hvec = -y;
                vvec = -z;
            }
            else {
                face = FACE_BACK;
                hvec = y;
                vvec = -z;
            }
        }
        else {
            if (z > 0) {
                face = FACE_TOP;
                hvec = -y;
                vvec = x;
            }
            else {
                face = FACE_BOTTOM;
                hvec = -y;
                vvec = -x;
            }
        }
    }
    else if (fabs(y) > fabs(z)) {
        if (y > 0) {
            face = FACE_LEFT;
            hvec = x;
            vvec = -z;
        }
        else {
            face = FACE_RIGHT;
            hvec = -x;
            vvec = -z;
        }
    }
    else {
        if (z > 0) {
            face = FACE_TOP;
            hvec = -y;
            vvec = x;
        }
        else {
            face = FACE_BOTTOM;
            hvec = -y;
            vvec = -x;
        }
    }

    image_x = (int)((face_dim / 2) + (face_dim / 2) * hvec);
    image_y = (int)((face_dim / 2) + (face_dim / 2) * vvec);

    /* ahem */
    if (image_x == face_dim)
        --image_x;
    if (image_y == face_dim)
        --image_y;

    assert(image_x >= 0);
    assert(image_x < face_dim);
    assert(image_y >= 0);
    assert(image_y < face_dim);

    /* Point to the pixel in which this star resides */
    pixel_dest = &faces[face][image_x + face_dim * image_y];

    /* Combine its brightness with any other star that
       happens to share the same pixel */
    *pixel_dest += (float) pixel_brightness(mag, mid_mag);
    if (*pixel_dest > 1) {
        *pixel_dest = 1;
        if (clipped != NULL)
            ++*clipped;
    }
}

int
alloc_faces(float **faces, int face_dim) {
    int f;

    for (f = 0; f < 6; ++f) {
        faces[f] = malloc(face_dim * face_dim * sizeof(float));
        if (faces[f] == NULL) {
            error(0, 0, "couldn't allocate memory for faces.\n");
            return -1;
        }
        memset(faces[f], 0, face_dim * face_dim * sizeof(float));
    }

    return 0;
}


/* `file' should be the file handle positioned at the start of the CSV file.
   `faces' should point to an array of six pointers to float, each of which
   should point to an array of face_dim * face_dim floats.
   ra_delta is added to to the right ascension of all the stars.
   Return the number of malformed lines in the CSV file.
   Stars with magnitude mid_mag are drawn as a middle grey (0.5). */
 
int
process_star_file(FILE *file, float **faces, int face_dim, double ra_delta,
        double dec_delta, double mid_mag, double cutoff_mag,
        int ra_field, int dec_field, int mag_field, int skip_lines) {
    int malformed_lines = 0;
    int c;
    int ret;
    int line_no = 1;
    int stars_processed = 0;
    int stars_drawn = 0;
    int stars_discarded = 0;
    int stars_processed_next_report = 100;
    int report_increment = 100;
    int removed_sun = 0;
    int clipped = 0;


    /* Read and discard the first N lines of the CSV file if necessary */
    while (skip_lines > 0) {
        while ((c = fgetc(file)) != '\n' && c != EOF);
        --skip_lines;
        ++line_no;
        if (c == EOF)
            break;
    }

    do {
        double ra, dec, mag;

        ret = read_star_line(file, &ra, &dec, &mag, line_no++,
                ra_field, dec_field, mag_field);

        if (ret == 1 && ra == 0 && dec == 0 && mag < -20) {
            /* If this is the Sun, don't draw it. */
            ++stars_processed;
            ++removed_sun;
        }
        else if (ret == 1 && mag > cutoff_mag) {
            /* Star is too dim; don't draw it. */
            ++stars_processed;
            ++stars_discarded;
        }
        else if (ret == 1) {
            draw_star(faces, ra, dec, mag, ra_delta, dec_delta,
                    face_dim, mid_mag, &clipped);

            ++stars_processed;
            ++stars_drawn;
        }
        else if (ret < 0) {
            ++malformed_lines;
        }

        /* Write a progress counter */
        if (stars_processed >= stars_processed_next_report) {
            progress_report(stars_processed, stars_drawn,
                    stars_discarded, malformed_lines,
                    removed_sun, clipped, 0);
            stars_processed_next_report += report_increment;
        }
    } while (ret != 0);

    progress_report(stars_processed, stars_drawn, stars_discarded,
            malformed_lines, removed_sun, clipped, 1);
    return malformed_lines;
}

int
output_faces(float **faces, int face_dim, int output_format) {
    int ret = 0;

    if (output_format == OFMT_PNM) {
        int f;
        for (f = 0; f < 6; ++f) {
            char filename[50];
            FILE *file;

            sprintf(filename, "face%d.pnm", f);

            file = fopen(filename, "w");
            if (file == NULL) {
                error(0, errno, "couldn't open %s for writing", filename);
            }
            else {
                float *face;
                int p;

                face = faces[f];

                fprintf(file, "P6\n");
                fprintf(file, "%d %d\n", face_dim, face_dim);
                fprintf(file, "%d\n", 255);

                for (p = 0; p < face_dim * face_dim; ++p) {
                    int val = face[p] * 255;
                    if (val > 255)
                        val = 255;
                    fputc(val, file);
                    fputc(val, file);
                    fputc(val, file);
                }

                if (fclose(file) == EOF) {
                    error(0, errno, "couldn't close %s", filename);
                    ret = 1;
                }
            }
        }
    }
    else {
        char *filename = "starbox.h";
        FILE *file = fopen(filename, "w");

        if (file == NULL) {
            error(0, errno, "Couldn't open starbox.h for writing");
        }
        else {
            int f;
            char *face_names[6] = { "front", "left", "back", "right", "top", "bottom"};
            /* Make header file idempotent */
            fprintf(file, "#ifndef _STARBOX_H\n");
            fprintf(file, "#define _STARBOX_H\n\n");
            for (f = 0; f < 6; ++f) {
                int r, c;

                fprintf(file, "float starbox_%s[%d][%d] = {\n",
                        face_names[f], face_dim,
                        face_dim);

                for (r = 0; r < face_dim; ++r) {
                    fprintf(file, "\t{");
                    for (c = 0; c < face_dim; ++c) {
                        if (c % 8 == 0)
                            fprintf(file, "\n\t\t");
                        fprintf(file, "%g, ", faces[f][r * face_dim + c]);
                    }
                    fprintf(file, "\n\t},\n");
                }
                fprintf(file, "};\n");
            }
            fprintf(file, "#endif\n");
            if (fclose(file) == EOF) {
                error(0, errno, "couldn't close %s", filename);
                ret = 1;
            }
        }
    }

    return ret;
}

float
get(float *arr, int dim, int x, int y) {
    // clip to image
    x = x<0 ? 0 : x;
    y = y<0 ? 0 : y;
    x = x>=dim ? dim-1 : x;
    y = y>=dim ? dim-1 : y;
    return arr[y*dim + x];
}

void
post_process_faces(float **faces, int face_dim) {
    for (int i=0 ; i<6 ; ++i) {
        float *src = faces[i];
        float *dst = malloc(face_dim * face_dim * sizeof(float));
        if (dst == NULL) {
            error(0, 0, "couldn't allocate memory for faces.\n");
            return;
        }
        // started off as a gaussian but was tweaked to look more appealing...
        float gaussian[5][5] = {
            {450, 26,  7,  0,  0 },
            {26, 16,  4,  0,  0 },
            { 7,  4,  1,  0,  0 },
            { 0,  0,  0,  0,  0 },
            { 0,  0,  0,  0,  0 },
        };
        float total = 0.0;
        for (int gy=-4 ; gy<=4 ; ++gy) {
            for (int gx=-4 ; gx<=4 ; ++gx) {
                int gy_ = gy>0 ? gy : -gy;
                int gx_ = gx>0 ? gx : -gx;
                total += gaussian[gy_][gx_];
            }
        }
        
        for (int y=0 ; y<face_dim ; ++y) {
            for (int x=0 ; x<face_dim ; ++x) {
                float v = 0;
                for (int gy=-4 ; gy<=4 ; ++gy) {
                    for (int gx=-4 ; gx<=4 ; ++gx) {
                        int gy_ = gy>0 ? gy : -gy;
                        int gx_ = gx>0 ? gx : -gx;
                        v += (gaussian[gy_][gx_]) * get(src, face_dim, x+gx, y+gy);
                    }
                }
                dst[y*face_dim + x] = v/total;
            }
        }
        free(faces[i]);
        faces[i] = dst;
    }
}

void
brighten_faces(float **faces, int face_dim, float scale_factor) {
    for (int i=0 ; i<6 ; ++i) {
        float *src = faces[i];
        for (int y=0 ; y<face_dim ; ++y) {
            for (int x=0 ; x<face_dim ; ++x) {
                src[y*face_dim + x] *= scale_factor;
            }
        }
    }
}

void
print_help(void) {
    printf("starbox - a program to take star position and magnitude data and produce the\n"
        "six faces of a skybox with the celestial sphere projected onto them.\n");
    printf("\n");
    printf("starbox [options] starfile.csv\n");
    printf("\t-h: show this help\n");
    printf("\t-v: version information\n");
    printf("\n");
    printf("\t-r <delta>: right ascension adjustment (degrees).\n");
    printf("\t-d <delta>: declination adjustment (degrees).\n");
    printf("\n");
    printf("\t-g <magnitude>: star magnitude to be drawn as middle grey (default %d).\n", DEFAULT_MID_MAG);
    printf("\t-C <cutoff>: don't draw stars with magnitude greater than <cutoff>.\n");
    printf("\n");
    printf("\t-p: output six PNM files called face<num>.pnm (default).\n");
    printf("\t    face0 is front, face1 is left, face2 is behind, face3 is\n"
        "\t    right, face4 is above (top row of image behind observer) and\n"
        "\t    face5 is below (top row of image in front of observer).\n");
    printf("\t-c: output C header file called starbox.h.\n");
    printf("\n");
    printf("\t-f <dimension>: height/width of a face in pixels (default %d).\n", DEFAULT_FACE_DIM);
    printf("\t-G: draw ascension and declination lines at 10 degree intervals.\n");
    printf("\n");
    printf("\t-1: skip first line of CSV file.\n");
    printf("\t-F <ra>,<dec>,<mag>: CSV field numbers for right ascension,\n"
        "\t    declination and magnitude (default 1,2,3).\n");
    printf("\t-H: equivalent to \"-F 8,9,14 -1\". Useful for reading HYG data.\n");
}

void
print_version(void) {
    printf("starbox - a program to produce a starfield texture for a skybox.\n");
    printf("Version 1.0, Graeme Cole 2010\n");
}

int main(int argc, char **argv) {
    FILE *starfile;
    int c;
    int ret;
    int malformed_lines = 0;
    float *faces[6];
    int output_format = OFMT_PNM;
    int face_dim = DEFAULT_FACE_DIM;
    double mid_mag = DEFAULT_MID_MAG;
    double ra_delta = 0, dec_delta = 0;
    float scale_factor = DEFAULT_SCALE_FACTOR;
    double cutoff_mag = DEFAULT_CUTOFF_MAG;
    int ra_field = 1, dec_field = 2, mag_field = 3;
    int skip_lines = 0;
    int draw_grid = 0;

    if (argc == 1) {
        print_help();
        exit(0);
    }

    while ((c = getopt(argc, argv, "1hHr:d:f:pcg:S:C:F:Gv")) != -1) {
        switch (c) {
            case '1':
                skip_lines++;
                break;
            case 'H':
                ra_field = 8;
                dec_field = 9;
                mag_field = 14;
                skip_lines = 1;
                break;
            case 'F':
                if (sscanf(optarg, "%d , %d , %d", &ra_field,
                        &dec_field, &mag_field) < 3) {
                    error(1, 0, "-F requires three comma-separated field numbers.");
                }
                if (ra_field < 1 || dec_field < 1 || mag_field < 1)
                    error(1, 0, "-F requires three positive field numbers.");
                break;
            case 'S':
                if (sscanf(optarg, "%f", &scale_factor) < 1) {
                    error(1, 0, "-S requires a numeric argument");
                }
                break;
            case 'C':
                if (sscanf(optarg, "%lf", &cutoff_mag) < 1) {
                    error(1, 0, "-C requires a numeric argument");
                }
                break;
            case 'f':
                if (sscanf(optarg, "%d", &face_dim) < 1 ||
                        face_dim <= 0) {
                    error(1, 0, "-f requires a positive argument");
                }
                break;
            case 'd':
                if (sscanf(optarg, "%lf", &dec_delta) < 1 || dec_delta < -180 || dec_delta > 180) {
                    error(1, 0, "-d requires a numeric argument between -180 and 180.");
                }
                break;
            case 'r':
                if (sscanf(optarg, "%lf", &ra_delta) < 1 || ra_delta < -360 || ra_delta > 360) {
                    error(1, 0, "-r requires a numeric argument between -360 and 360");
                }
                break;
            case 'h':
                print_help();
                exit(0);

            case 'v':
                print_version();
                exit(0);
            
            case 'p':
                output_format = OFMT_PNM;
                break;
            case 'c':
                output_format = OFMT_C;
                break;

            case 'g':
                if (sscanf(optarg, "%lf", &mid_mag) < 1) {
                    error(1, 0, "-g requires a numeric argument.");
                }
                break;

            case 'G':
                draw_grid = 1;
                break;

            default:
                error(1, 0, "Unknown option -%c\n", optopt);
        }
    }

    /* Optionless argument is the star file */
    if (optind >= argc) {
        starfile = stdin;
    }
    else {
        char *filename = argv[optind];
        if (!strcmp(filename, "-")) {
            starfile = stdin;
        }
        else {
            starfile = fopen(filename, "r");
            if (starfile == NULL) {
                error(1, errno, "Couldn't open %s\n", filename);
            }
        }

        /* Check for any junk options */
        if (optind + 1 < argc) {
            error(1, 0, "Unexpected argument %s\n", argv[optind + 1]);
        }
    }

    /* Allocate and intialise the faces */
    if (alloc_faces(faces, face_dim) == -1)
        exit(1);

    /* Draw grid if necessary */
    if (draw_grid) {
        double ra, dec, mag;

        fprintf(stderr, "Drawing grid...");
        for (dec = -89; dec <= 89; ++dec) {
            for (ra = 0; ra < 360; ++ra) {
                if (dec == 0 || ra == 0)
                    mag = 3.3;
                else
                    mag = 4;
                if ((int) ra % 10 == 0 || (int) abs(dec) % 10 == 0) {
                    draw_star(faces, ra, dec, mag, ra_delta, dec_delta,
                        face_dim, 2.5, NULL);
                }
            }
        }
        fprintf(stderr, " done.\n");
    }

    /* Process all the stars in the file */
    malformed_lines = process_star_file(starfile, faces, face_dim, ra_delta,
            dec_delta, mid_mag, cutoff_mag, ra_field, dec_field,
            mag_field, skip_lines);

    fclose(starfile);

    if (malformed_lines > 0) {
        fprintf(stderr, "%d malformed lines.\n", malformed_lines);
    }

    fprintf(stderr, "Post processing...\n");
    post_process_faces(faces, face_dim);
    brighten_faces(faces, face_dim, scale_factor);

    /* Output the results in whatever format the user wanted */
    ret = output_faces(faces, face_dim, output_format);

    /* Clean up */
    for (c = 0; c < 6; ++c)
        free(faces[c]);

    /* Clear off */
    if (ret == 0)
        return 0;
    else
        return 1;
}


