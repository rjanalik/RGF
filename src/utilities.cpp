#include "utilities.hpp"
#include <algorithm>
#include <iomanip>
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
struct option {
    option(std::string short_name, std::string long_name, std::string type, std::string description)
        : short_name_(short_name), long_name_(long_name), type_(type), description_(description){};
    std::string short_name_;
    std::string long_name_;
    std::string type_;
    std::string description_;
    void print_message(std::ostream &stream = std::cout) {
        stream << "-" + short_name_ + "," << std::setw(20) << "-" + long_name_ << std::setw(50) << "[" + type_ + "]" + "\n"
               << std::setw(15) << description_ << std::endl;
    }
};

void print_help_message(const std::vector<option> &options, std::ostream &stream = std::cout,
                        std::string const &which = std::string()) {
    if (!which.empty()) {

    } else {
        for (auto opt : options) {
            opt.print_message();
        }
    }
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) { return std::find(begin, end, option) != end; }

char *getCmdOption(char **begin, char **end, const std::string &option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

void parse_args(int argc, char *argv[], std::string &base_path, size_t &ns, size_t &nt, size_t &nb, size_t &no) {
    std::vector<option> options = {
        option("-h", "--no", "std::string", "print this help message"),
        option("-p", "--path", "std::string", "base path to folder containing all files"),
        option("-s", "--ns", "integer", "number of spatial effects"),
        option("-t", "--nt", "integer", "number of temporal effects"),
        option("-f", "--nb", "integer", "number of fixed effects"),
        option("-n", "--no", "integer", "number of data samples"),
    };
    if (cmdOptionExists(argv, argv + argc, "-h") || cmdOptionExists(argv, argv + argc, "--help")) {
        print_help_message(options);
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
