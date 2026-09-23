/* This file is part of the AENET package.
 *
 * Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
 *
 * This Source Code Form is subject to the terms of the Mozilla Public License,
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can
 * obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "aenet.h"
#include <dlfcn.h>
#include <stdio.h>

typedef void (*init_function)(int, char *[], int *);
typedef void (*final_function)(int *);
typedef AENET_BOOL (*loaded_function)(void);
typedef void (*convert_function)(int, char *[], int, int [], int [], int *);

static void *symbol(void *library, const char *name)
{
    void *value = dlsym(library, name);
    if (value == NULL) {
        fprintf(stderr, "missing C API symbol %s: %s\n", name, dlerror());
    }
    return value;
}

int main(int argc, char **argv)
{
    void *library;
    init_function init;
    final_function final;
    loaded_function all_loaded;
    convert_function convert;
    int *ok;
    char *types[] = {"Cu", "Au"};
    char *incoming[] = {"Au", "Cu"};
    int input[] = {1, 2};
    int output[] = {0, 0};
    int status = -1;

    if (argc != 2) {
        fprintf(stderr, "usage: %s /path/to/libaenet.so\n", argv[0]);
        return 64;
    }
    library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    if (library == NULL) {
        fprintf(stderr, "cannot load %s: %s\n", argv[1], dlerror());
        return 1;
    }
    init = (init_function)symbol(library, "aenet_init");
    final = (final_function)symbol(library, "aenet_final");
    all_loaded = (loaded_function)symbol(library, "aenet_all_loaded");
    convert = (convert_function)symbol(library, "aenet_convert_atom_types");
    ok = (int *)symbol(library, "AENET_OK");
    if (init == NULL || final == NULL || all_loaded == NULL ||
        convert == NULL || ok == NULL) {
        return 2;
    }
    init(2, types, &status);
    if (status != *ok || all_loaded()) {
        return 3;
    }
    convert(2, incoming, 2, input, output, &status);
    if (status != *ok || output[0] != 2 || output[1] != 1) {
        return 4;
    }
    final(&status);
    if (status != *ok) {
        return 5;
    }
    dlclose(library);
    puts("C API initialization, type mapping, and finalization passed");
    return 0;
}
