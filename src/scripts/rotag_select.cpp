#include <getopt.h>
#include <iostream>

int main(int argc, char *argv[]) {
  const struct option longopts[] = {
    {"target",         0, 0, 't'},
    {"select",         0, 0, 's'},
    {"tags",           0, 0     },
    {"related-data",   0, 0, 'r'},
    {"pdb",            0, 0, 'p'},
    {"keep-ignored",   0, 0, 'k'},
    {"random-seed",    0, 0, 'x'},
    {"help",           0, 0, 'h'},
    {"version",        0, 0, 'v'},
  };

  int index;
  int iarg = 0;

  while(iarg != -1) {
    iarg = getopt_long(argc, argv, "s:vh", longopts, &index);
  }

  return 1;
}
