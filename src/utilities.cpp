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
template <class option_type> struct option {
  option(std::string short_name, std::string long_name, std::string type,
         std::string description, bool required)
      : short_name_(short_name), long_name_(long_name), type_(type),
        description_(description), is_required_(required){};
  option(std::string short_name, std::string long_name, std::string type,
         std::string description, bool required, option_type *value)
      : short_name_(short_name), long_name_(long_name), type_(type),
        description_(description), is_required_(required) {
    value_ptr_ = value;
  };
  std::string short_name_;
  std::string long_name_;
  std::string type_;
  std::string description_;
  option_type *value_ptr_ = nullptr;
  bool is_required_;
  void print_message(std::ostream &stream = std::cout) const {
    stream << short_name_ + ", " + long_name_ << std::right << std::setw(60)
           << "[" + type_ + "]" << std::endl
           << "\t" << description_ << "\n\n";
  }
};

template <class T>
void print_help_message(const option<T> &option,
                        std::ostream &stream = std::cout) {
  option.print_message();
}

template <class T>
void print_help_message(const std::vector<option<T>> &options,
                        std::ostream &stream = std::cout) {
  for (auto opt : options) {
    opt.print_message();
  }
}

template <class T>
bool is_option_present(char *begin[], char *end[], T const &option) {
  return (is_option_present(begin, end, option.short_name_) ||
          is_option_present(begin, end, option.long_name_));
}

template <>
bool is_option_present(char *begin[], char *end[], const std::string &option) {
  return std::find(begin, end, option) != end;
}

template <class T>
void set_cml_option(char *begin[], char *end[], option<T> &option) {
  char **itr;
  itr = std::find(begin, end, option.short_name_);
  if (itr != end && ++itr != end)
    *(option.value_ptr_) = *itr;
  else {
    itr = std::find(begin, end, option.long_name_);
    if (itr != end && ++itr != end)
      *(option.value_ptr_) = *itr;
  }
}
template <>
void set_cml_option(char *begin[], char *end[], option<size_t> &option) {
  char **itr;
  itr = std::find(begin, end, option.short_name_);
  if (itr != end && ++itr != end)
    *(option.value_ptr_) = atoi(*itr);
  else {
    itr = std::find(begin, end, option.long_name_);
    if (itr != end && ++itr != end)
      *(option.value_ptr_) = atoi(*itr);
  }
}

template <class T>
void check_for_missing_option(char **begin, char **end,
                              std::vector<option<T>> const &options,
                              std::ostream &stream = std::cout) {
  for (auto opt : options) {
    if (opt.is_required_)
      if (!is_option_present(begin, end, opt)) {
        stream << "Missing Option:\n";
        print_help_message(opt);
      }
  }
}

void parse_args(int argc, char *argv[], std::string &base_path, size_t &ns,
                size_t &nt, size_t &nb, size_t &no, std::ostream &stream) {
  struct option<std::string> help_option = {
    "-h", "--help", "std::string", "print this help message", false
  };
  std::vector<option<std::string>> string_options = {
      help_option,
      {"-p", "--path", "std::string",
       "base path to folder containing all files", true, &base_path}};
  std::vector<option<size_t>> number_options = {
      {"-s", "--ns", "integer", "number of spatial effects", true, &ns},
      {"-t", "--nt", "integer", "number of temporal effects", true, &nt},
      {"-f", "--nb", "integer", "number of fixed effects", true, &nb},
      {"-n", "--no", "integer", "number of data samples", true, &no},
  };
  if (is_option_present(argv, argv + argc, help_option)) {
    print_help_message(string_options);
    print_help_message(number_options);
    exit(0);
  }
  check_for_missing_option(argv, argv + argc, string_options);
  check_for_missing_option(argv, argv + argc, number_options);
  for (auto opt : string_options) {
    set_cml_option<std::string>(argv, argv + argc, opt);
  }
  for (auto opt : number_options) {
    set_cml_option<size_t>(argv, argv + argc, opt);
  }
}

namespace utilities {
  void print_header(std::string title, size_t length, char symbol, size_t sep_width, char sep, std::ostream &stream) {
      size_t title_length = title.length();
      size_t symbol_length = length-title.length()-2*sep_width;
      // print only header line
      if(title_length == 0){
          stream << std::string(length, symbol) << std::endl;
      } else if(symbol_length <= 0){
          // print only title if title to long
          stream << title << std::endl;
      } else{
          stream << std::string(symbol_length/2, symbol) << std::string(sep_width, sep) << title << std::string(sep_width, sep) << std::string(symbol_length/2, symbol) << std::endl;
      }
  }
}
