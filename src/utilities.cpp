#include "utilities.hpp"
#include <algorithm>
#include <iostream>
#include <string>

void icopy(int N, int *x, int *x_copy) {
  int i;

  for (i = 0; i < N; i++) {
    x_copy[i] = x[i];
  }
}

/************************************************************************************************/

int Round(double x) {
  if (x < 0) {
    return (int)(x - 0.5);
  } else {
    return (int)(x + 0.5);
  }
}

/************************************************************************************************/

double get_time(double t0) {
  timeval tim;
  gettimeofday(&tim, NULL);

  return (double)(tim.tv_sec + (tim.tv_usec / 1000000.0)) - t0;
}

/************************************************************************************************/
void print_help_message(std::ostream &stream = std::cout) {
  stream << "RGF Call : path_to_folder ns nt nb no" << std::endl;

  stream << "[string:base_path]          path to folder containing "
            "all files "
         << std::endl;
  stream << "[integer:ns]                number of spatial "
            "grid points "
         << std::endl;
  stream << "[integer:nt]                number of temporal "
            "grid points "
         << std::endl;
  stream << "[integer:nb]                number of fixed effects" << std::endl;

  stream << "[integer:no]                number of data samples" << std::endl;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}
void parse_args(int argc, char *argv[], std::string &base_path, size_t &ns,
                size_t &nt, size_t &nb, size_t &no) {
  if (cmdOptionExists(argv, argv + argc, "-h") ||
      cmdOptionExists(argv, argv + argc, "--help")) {
    print_help_message();
    exit(1);
  }
  if (cmdOptionExists(argv, argv + argc, "-p")) {
    base_path = getCmdOption(argv, argv + argc, "-p");
  } else if (cmdOptionExists(argv, argv + argc, "--path")) {
    base_path = getCmdOption(argv, argv + argc, "--path");
  }
  if (cmdOptionExists(argv, argv + argc, "-s")) {
    ns = atoi(getCmdOption(argv, argv + argc, "-s"));
  } else if (cmdOptionExists(argv, argv + argc, "--ns")) {
    ns = atoi(getCmdOption(argv, argv + argc, "--ns"));
  }
  if (cmdOptionExists(argv, argv + argc, "-t")) {
    nt = atoi(getCmdOption(argv, argv + argc, "-t"));
  } else if (cmdOptionExists(argv, argv + argc, "--nt")) {
    nt = atoi(getCmdOption(argv, argv + argc, "--nt"));
  }
  if (cmdOptionExists(argv, argv + argc, "-f")) {
    nb = atoi(getCmdOption(argv, argv + argc, "-f"));
  } else if (cmdOptionExists(argv, argv + argc, "--nb")) {
    nb = atoi(getCmdOption(argv, argv + argc, "--nb"));
  }
  if (cmdOptionExists(argv, argv + argc, "-n")) {
    no = atoi(getCmdOption(argv, argv + argc, "-n"));
  } else if (cmdOptionExists(argv, argv + argc, "--no")) {
    no = atoi(getCmdOption(argv, argv + argc, "--no"));
  }
}
