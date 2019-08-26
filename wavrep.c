// wavrep.c - copy WAV file with repetition

#include <stdio.h>
#include <stdlib.h>
#include <sigpro.h>

static char *ifn, *ofn;
static int nreps = 1;
static int nbits = 16;

static void
usage()
{
    fprintf(stdout, "usage: wavrep [options] <input_file> <output_file>\n");
    fprintf(stdout, "options\n");
    fprintf(stdout, "-b N  set number of bits to N\n");
    fprintf(stdout, "-n N  set number of repetitions to N\n");
    exit(0);
}

static void
parse_args(int ac, char *av[])
{
    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'b' && ac > 2) {
                nbits = atoi(av[2]);
                ac--;
                av++;
            } else if (av[1][1] == 'n' && ac > 2) {
                nreps = atoi(av[2]);
                ac--;
                av++;
            }
            ac--;
            av++;
        } else {
            break;
        }
    }
    if (ac < 3) usage();
    ifn = av[1];
    ofn = av[2];
}

int
main(int ac, char *av[])
{
    float fs;
    static VAR *vl;

    parse_args(ac, av);
    vl = sp_wav_read(ifn, 0, 0, &fs);
    fprintf(stdout, " input file: %s fs=%.0f\n", ifn, fs);
    fprintf(stdout, "output file: %s reps=%d bits=%d\n", ofn, nreps, nbits);
    sp_wav_write_rep(ofn, vl, &fs, nbits, nreps);
    sp_var_clear(vl);
}
