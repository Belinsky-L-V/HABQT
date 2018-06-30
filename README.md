# HABQT

Aim of the project is to extend adaptive Bayesian quantum tomography as
described in [2012 paper by Huszár and
Houlsby](https://doi.org/10.1103/PhysRevA.85.052120) by using a hierarchical
distribution over density matrices of all possible ranks.

Includes:

* a Haskell library
* a shared library which provides a C interface to the tomography function
* an executable that simulates tomography of random states and outputs
  infidelity between true states and mean Bayesian estimates to a file

Please refer to `HABQT-simulation --help` for executable usage instructions,
accompanying Haddock documentation for Haskell API, and [libHABQT header
file](./libHABQT.h) for shared library C ABI description.

### Installation instructions

#### Linux

No special setup should be necessary. Simply make sure GSL, BLAS and LAPACK are
installed on your system and install normally using stack, cabal or Setup.hs.

#### Windows

Making the necessary shared libraries and tools available on windows can be a
bit tricky. One way to do this is install them inside MSYS2 that comes with
stack on windows (on x86\_64 can be found under
`%LOCALAPPDATA%\Programs\stack\x86_64-windows`):

0. MSYS may not be present, in which case execute `stack build` in HABQT
   directory. The build will fail due to missing libraries/tools, but MSYS
   should be installed after it, and you will be able to add them.
1. Launch msys2.exe
2. (Optional) Update MSYS2 with `pacman -Syu`. It may be necessary to restart
   the shell: follow the instructions displayed in it.
3. Install the appropriate mingw toolchain, which includes necessary tools like
   pkg-config. E.g.  
   ```pacman -S mingw-w64-x86_64-toolchain```
4. Install GSL for the appropriate mingw toolchain (the 64-bit one this case):  
   ```pacman -S mingw64/mingw-w64-x86_64-gsl```
5. Install openblas for appropriate mingw toolchain:  
   ```pacman -S mingw64/mingw-w64-x86_64-openblas```
6. The versions/naming conventions hmatrix expects differ from what is used in
   modern MSYS2, so it's necessary to either link of create renamed copies of
   two libraries from the appropriate mingw toolchain (in my case found under
   `%LOCALAPPDATA%\Programs\stack\x86_64-windows\msys2-20150512\mingw64\bin`):
   `libgfortran-4.dll` to `libgfortran.dll` and `libgslcblas-0.dll` to
   `gsl-0.dll`. I recommend placing links/copies in some directory that isn't
   normally on PATH (neither windows nor MSYS) and explicitly pointing stack to
   them during installation with `--extra-lib-dirs ` (example follows in next
   step).
7. Outside MSYS2 open a normal windows shell, navigate to HABQT folder and
   build/install with stack, passing appropriate flags and library dirs:  
   ```stack build --flag hmatrix:openblas --extra-lib-dirs=D:\lib```  
   where `D:\lib` contains ` libgfortran.dll` and `gsl-0.dll`.

To use the shared library or executable, you’ll need to have several mingw
libraries on PATH or in the same directory: libgcc_s_seh-1.dll,
libgfortran-4.dll, libgsl-23.dll, libgslcblas-0.dll, libopenblas.dll,
libquadmath-0.dll and libwinpthread-1.dll.
