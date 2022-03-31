# RGF
    .
    ├── docs                    # Documentation with `index.html`
    ├── build                   # Build files
    │   ├── obj                 # Object Files
    │   └── bin                 # Binary 
    ├── data                    # Data files 
    │   ├── input               # Input files/matrices
    │   │   │── ghcn            # Global Historical Climatology Network (GHCN)
    │   │   │    ├── 2018 
    │   │   │    ├── 2019
    │   │   │    └── ...
    │   │   │── Radim           # Original test files 
    │   │   │
    │   │   │
    │   │   └── tests           # test files created with helper script
    │   │
    │   └── output              # Output Files From ./main
    ├── applications            # External Libraries s.a. magma
    ├── include                 # Header Files
    ├── scripts                 # Scripts i.e. for setting up the environment
    │   ├── R                   # R script, mainly for data generation of ghcn
    │   ├── C                   # C script, for data generation of test data
    │   ├── setup_kaust.sh      # C script, for data generation of test data
    │   ├── C                   # C script, for data generation of test data
    │   └── run_single.sh
    └── src                    # Source files
## TODO:
- Check version conflicts:  
`warning: libcusparse.so.11, needed by /home/x_pollakgr/applications/magma-2.5.4/lib/libmagma.so, may conflict with libcusparse.so.10`
- CUDA header dependencies not in makefile included
## Installation
## Fetching the data & Creating the Model
### Fetching the data 
1. Change to:
    ``` sh
    <project root>/scripts/R
    ```
2. Run: 
   ```sh
   Rscript fetch_ghcn_daily.R [--help]
   ```
### Creating The Models
1. Change to:
    ``` sh
    <project root>/scripts/R
    ```
2. Run: 
   ```sh
   Rscript generate_spatial_temporal_model.R [--help]
   ```
### Compiling The Code
#### KAUST
1. Setup the enviroment:
- Run the setup script:
   ```shell
   source ~/RGF/scripts/setup_kaust.sh
   ```
2. Run the makefile:
    ``` sh
    make RGF_VERSION=<BASE,ASYNC,ASYNC_2S,BANDED,PARDISO> [DEBUG=Level i.e. 1 or 2] -B 
    ```
### Running the code 
#### KAUST
##### To run single test file: 
Helper script to run the code, hard coded for test files (adapt as fit):
    ```sh
    ./scripts/run_single.sh ns nt nb
    ```
##### To run single test file and profile it using nvvp/nsys (comment/uncomment line 34)
Helper script to profile the code, hard coded for test files (adapt as fit):
    ```sh
    sudo ./scripts/run_single.sh ns nt nb nvvp_file_name
    ```
##### To run all test files from a folder using the previous script:
    to run all version from a given folder:
    ```sh
    python ./scripts/run_all.py
    ```
### Generating Examples
    From the folder `RGF/scripts/C` run
    ```sh
    bash ./generate_tests.sh ns nt nb
    ```
### Generating and Running
    From the folder `RGF/scripts/C` run
    ```sh
    bash ./generate_and_run ns nt nb
    ```
