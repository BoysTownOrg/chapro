// wavrep.c - copy WAV file with repetition

#include <stdio.h>
#include <stdlib.h>
#include <sigpro.h>

static void
usage()
{
    fprintf(stdout, "usage: wavrep [-n <count>] input_file output_file\n");
    exit(0);
}

int
main(int ac, char *av[])
{
    char *ifn, *ofn;
    float fs;
    int nreps = 1;
    int nbits = 16;
    static VAR *vl;

    while (ac > 1) {
        if (av[1][0] == '-') {
            if (av[1][1] == 'n' && ac > 2) {
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
    vl = sp_wav_read(ifn, 0, 0, &fs);
    fprintf(stdout, " input file: %s fs=%.0f\n", ifn, fs);
    fprintf(stdout, "output file: %s reps=%d\n", ofn, nreps);
    sp_wav_write_rep(ofn, vl, &fs, nbits, nreps);
    sp_var_clear(vl);
}
