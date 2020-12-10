#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#define PUBLIC_VARS 3
#define SECRET_VARS 3

static int KEY = 7;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Wrong number of arguments\n");
        return 1;
    }

    int input = atoi(argv[1]);

    bool v[PUBLIC_VARS] = { 0 };
    for (int i = 0; i < PUBLIC_VARS; ++i) {
        v[i] = input & (0x1 << i);
    }

    // preprocessing phase
    if (argc == 3) {
        KEY = atoi(argv[2]);
    }
    bool x[SECRET_VARS] = { 0 };
    for (int i = 0; i < SECRET_VARS; ++i) {
        x[i] = KEY & (0x1 << i);
    }

    bool bool_function = (v[0] && v[1] && x[0]) ^ (v[0] && v[1] && x[1]) ^ (v[2] && x[0] && x[2]) ^ (v[1] && x[2]) ^ (v[0] && x[0]) ^ (v[0] && v[1]) ^ (x[0] && x[2]) ^ (v[1]) ^ (x[2]) ^ (true);

    return bool_function;
}
