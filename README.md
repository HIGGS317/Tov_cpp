# Tov_cpp
C++ solver for Tolman–Oppenheimer–Volkoff (TOV)  Equation for Neutron Stars

This solver generates and writes mass and radius sequence in a csv-like file for a given equation of state.
The solver can also calculate love-number and compactness based on user's input
at runtime.

``` 
1. Clone the repo
2. cd to the cloned directory
```

# Build 

```cmake
cmake .
make
```
This will create a **bin** folder in the current directory containing
executable named Tov_cpp.
# Run

````shell
cd bin
./Tov_cpp [path to eos file] [number of points]
````
# Important points to note
1. Make sure that eos file is a csv file
2. Eos file should contain at least two columns of pressure and energy density
3. Path to eos file should be relative to executable.
