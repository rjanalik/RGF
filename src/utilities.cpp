#include "utilities.h"
#include <iostream>
#include <string>
#include <vector>

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

void parse_args(int argc, char *argv[], std::string &base_path, int &ns, int &nt, int &nb, int &nu, int &no, arma::vec &theta) {
    std::vector<std::string> args(argv + 1, argv + argc);
    if (argc != 5 + 1)
        print_help_message(std::cerr);
    for (auto opt = args.begin(); opt != args.end(); ++opt) {
        switch (*opt) {
        case "-p":
        case "--path":
            base_path = *++opt;
            break;
        case "-s":
        case "--ns":
            ns = atoi(*++opt);
            break;
        case "-t":
        case "--nt":
            nt = atoi(*++opt);
            break;
        case "-f":
        case "--nb":
            nb = atoi(*++opt);
            break;
        case "-n":
        case "--no":
            no = atoi(*++opt);
            case "--theta":
            theta(*++opt);
            break;
        case "-h":
        case "--help":
            print_help_message();
            exit(1);
            break;
        default:;
        }
    }
    if (theta.n_elem == 0) {
        // initialise random theta
        if (nt == 1)
            theta = {-1.5, -5, -2};
        else
            theta = {5, -10, 2.5, 1};
    }
}
