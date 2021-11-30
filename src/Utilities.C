#include "Utilities.H"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>


#if 0
typedef CPX T;
#define assign_T(val) CPX(val, 0.0)
#else
typedef double T;
#define assign_T(val) val
#endif

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
                size_t &nt, size_t &nb, std::ostream &stream) {
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

  void read_test_matrix_nnz(size_t &nnz, std::string file_path){
    FILE *F = fopen(file_path.c_str(),"r");
    fscanf(F,"%zu",&nnz);
    fscanf(F,"%zu",&nnz);
    fscanf(F,"%zu",&nnz);
    fclose(F);
  }
  void read_test_matrix(size_t *ia, size_t *ja, double *a, double *rhs, size_t rows, size_t nnz, size_t nrhs, std::string mat_path, std::string rhs_path){
    FILE *F = fopen(mat_path.c_str(), "r");
    double val;
    size_t f_rows;
    size_t f_nnz;
    int i;
    /* read in matrix A, sparse matrix in CSR format */
    fscanf(F,"%zu",&f_rows);
    fscanf(F,"%zu",&f_rows);
    fscanf(F,"%zu",&f_nnz);
    for (i = 0; i <= rows; i++){
       fscanf(F,"%zu",&ia[i]);
    }
    for (i = 0; i < ia[rows]; i++){
       fscanf(F,"%zu",&ja[i]);
    }
    for (i = 0; i < ia[rows]; i++){
       fscanf(F,"%lf",&val);
       a[i] = assign_T(val);
    }
    fclose(F);
    // READ RHS
    nrhs   = 1;
    rhs    = new T[nrhs*rows];

    F = fopen(rhs_path.c_str(), "r");
    for (int i = 0; i < nrhs*rows; i++)
    {
        fscanf(F,"%lf",&rhs[i]);
    }
    fclose(F);
  }

  void print_header(std::string title, size_t length, char symbol,
                    size_t sep_width, char sep, std::ostream &stream) {
    size_t title_length = title.length();
    size_t symbol_length = length - title.length() - 2 * sep_width;
    // print only header line
    if (title_length == 0) {
      stream << std::string(length, symbol) << std::endl;
    } else if (symbol_length <= 0) {
      // print only title if title to long
      stream << title << std::endl;
    } else {
      stream << std::string(symbol_length / 2, symbol)
            << std::string(sep_width, sep) << title
            << std::string(sep_width, sep)
            << std::string(symbol_length / 2, symbol) << std::endl;
    }
  }
//   enum Code {
//   BG_BLACK    = 40,
//   BG_RED      = 41,
//   BG_GREEN    = 42,
//   BG_YELLOW   = 43,
//   BG_BLUE     = 44,
//   BG_MAGENT   = 45,
//   BG_CYAN     = 46,
//   BG_WHITE    = 47,
//   BG_DEFAULT  = 49,
//   FG_BLACK    = 30,
//   FG_RED      = 31,
//   FG_GREEN    = 32,
//   FG_YELLOW   = 33,
//   FG_BLUE     = 34,
//   FG_MAGENT   = 35,
//   FG_CYAN     = 36,
//   FG_WHITE    = 37,
//   FG_DEFAULT  = 39
// };
//   class Modifier {
//     Code code;
//     public:
//       Modifier(Code pCode) : code(pCode) {}
//       friend std::ostream&
//       operator<<(std::ostream& os, const Modifier& mod) {
//         return os << "\033[" << mod.code << "m";
//       }
//   };
// }

  constexpr size_t nmax {200};

  size_t number_of_digits(double n) {
    std::ostringstream strs;
    strs << n;
    size_t digits = strs.str().size();
    return n < 0 ? digits+1 : digits;
  }

  void get_max_digits_per_column(const T *M, size_t n, size_t *max_len_per_col){
    size_t max_len;
    size_t num_len;
    for (size_t j = 0; j < n; ++j) {
      max_len = 0;
      for (size_t i = 0; i < n; ++i){
        num_len = number_of_digits(M[i*n+j]);
        if (max_len < num_len)
          max_len = num_len ;
      }
      max_len_per_col[j] = max_len+1;
    }
  }


  void print_matrix(T *M, size_t m, size_t n, bool RMO) {
    size_t *max_len_per_col = new size_t[n];
    get_max_digits_per_column(M, n, max_len_per_col);
    size_t max_len = *std::max_element(max_len_per_col , max_len_per_col + n);
    size_t idx;
    for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < n; ++j){
        if(RMO)
          idx = i*m+j;
        else
          idx = i+m*j;
        std::cout << (j == 0 ? "\n| " : "") << std::setw(max_len) << M[idx] << (j == n - 1 ? " |" : " ");
      }
        // colors \033[0;41m
    delete[] max_len_per_col;
    std::cout << '\n';
  }

  void print_matrix_structure(T *M, size_t m, size_t n, bool RMO) {
    std::string s;
    size_t idx;
    for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < n; ++j){
        if(RMO)
          idx = i*m+j;
        else
          idx = i+m*j;
        if (M[idx] < 0) {
          s = "-";
        } else if (M[idx] > 0){
          s = "+";
        } else{
          s = "0";
        }
        std::cout << (j == 0 ? "\n| " : "") << std::setw(3) << s << (j == n - 1 ? " |" : " ");
        // colors \033[0;41m
      }
    std::cout << '\n';
  }

  void print_matrix_structure_with_fixed_effects(T *M, size_t m, size_t n, size_t nd, bool RMO) {
    std::string s;
    size_t idx;
    for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < n; ++j){
        if(RMO)
          idx = i*m+j;
        else
          idx = i+m*j;
        if (j > n-nd-1 || i > n-nd-1) {
          s = "=";
        } else if (M[idx] < 0) {
          s = "-";
        } else if (M[idx] > 0){
          s = "+";
        } else{
          s = "0";
        }
        std::cout << (j == 0 ? "\n| " : "") << std::setw(3) << s << (j == n - 1 ? " |" : " ");
        // colors \033[0;41m
      }
    std::cout << '\n';
  }


  size_t mf_block_index(size_t supernode, size_t *diag_pos, size_t b_size) {
    return diag_pos[supernode * b_size];
  }

  size_t mf_block_lda(size_t supernode, size_t ns, size_t nt, size_t nd) {
    // return matrix->index_i[c*b_size];
    // two blocks
    if (supernode < nt - 1)
      return 2 * ns + nd;
    // one block
    if (supernode < nt)
      return ns + nd;
    // dense block
    else
      return nd;
  }

  void print_MF_structure(T *MF, size_t ns, size_t nt, size_t nd, size_t *diag_pos) {
    std::string s;
    size_t supenode_bs = 2*ns*ns;
    size_t fixed_effects_bs = nd*ns;
    size_t total_rows_of_MF = nt*ns+nd;
    size_t total_cols_of_MF = nt*ns+nd;
    size_t col_size = ns;
    // Initalize memory to zero
    T *M_out = new T[total_rows_of_MF*total_cols_of_MF]();
    std::cout << "nt-1 = " << nt-1 << std::endl;
    std::cout << "total_rows_of_MF = " << total_rows_of_MF << std::endl;
    std::cout << "col_size = " << col_size << std::endl;
    for (size_t supernode = 0; supernode < nt-1; supernode++) {
      size_t supernode_fc = supernode * ns;
      size_t supernode_lc = supernode < nt
                                ? (supernode + 1) * ns
                                : ns * nt + nd;
      size_t rows = mf_block_lda(supernode, ns, nt, nd);
      size_t ind = mf_block_index(supernode, diag_pos, ns);
      size_t cols = supernode_lc - supernode_fc;
      memcpy(&M_out[supernode*(nt*ns+nd)*ns + supernode*ns*ns], &MF[ind], rows*cols*sizeof(T));
    }
    print_matrix_structure(M_out, total_rows_of_MF, total_cols_of_MF);
    // for (size_t IB = 0; IB < nt - 1; IB++) {
    //   // size_t diag_block_start = diag_pos[IB * ns];
    //     for (size_t i = 0; i < total_rows_of_MF; i++) {
    //       for (size_t j = 0; j < col_size; j++) {
    //         // glb_idx = IB*(total_rows_of_MF*col_size);
    //         // if(IB*ns < i < IB*ns+2*ns){
    //         //   MF[diag_pos[IB*ns]];
    //         // }
    //       }
    //     }
    // }
  }


  void print_csr(size_t *ia, size_t *ja, double *a, size_t n, size_t nd, bool structureOnly){
    T *M = new T[n*n]();
    size_t elem = 0;
    for (size_t j = 0; j < n; ++j){
      for (size_t i = 0; i < (ia[j+1]-ia[j]); ++i){
        size_t row_idx = ja[ia[j]+i];
        // a stored as column major
        T val = a[elem++];
        M[row_idx+n*j] = val;
        if(structureOnly)
          M[row_idx*n+j] = val;
      }
    }
    if(structureOnly){
      if(nd>0)
        print_matrix_structure_with_fixed_effects(M, n, n, nd, false);
      else
        print_matrix_structure(M, n, n, false);
    }
    else
      print_matrix(M, n, n, false);
    delete[] M;
  }

  void print_ia_ja_a(const size_t *ia, const size_t *ja, const T *a, size_t n){
    size_t nnz = ia[n];
    std::cout << "ia:" << std::setw(10) << " ";
    std::cout << "ja:" << std::setw(10) << " ";
    std::cout << "a:" << std::setw(10) << "\n";
    for(size_t i = 0; i < nnz; i++){
      if(i < n){
        std::cout << ia[i] << std::setw(10) << " ";
      } else
        std::cout << std::setw(10) << " ";

      std::cout << ja[i] << std::setw(10) << " ";
      std::cout << a[i] << std::setw(10) << " ";
      std::cout << "\n";
    }
  }

  const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
  }
  bool file_exists(const std::string &file_name) { return std::fstream{file_name} ? true : false; }
  void if_not_exists_abort(std::string const file_name) {
      if(file_exists(file_name))
          return;
      std::cerr << file_name
                << " couldn\'t be opened (not existing or failed to open)\n";
      exit(1);
  }
  void if_not_exists_abort(std::initializer_list<std::string> const file_names) {
      for (std::string file_name : file_names)
          if_not_exists_abort(file_name);
  }

} // namespace utilities
