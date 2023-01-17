# OpenOptions

This package, "OpenOptions", is an R package that provides a standard set of functions for evaluating options according to the Binomial and Black-Sholes-Merton models. The package is written in C++ for efficient implementation of option pricing, calculating Greeks, and implied volatilities.

This package is useful for finance professionals, researchers, and students who are interested in option pricing and related calculations. The package is easy to use and provides efficient and accurate results.

# Package Installtion 

## Installation via R

1. Open R on your computer.

2. Install the "devtools" package by running the command:

```{r}
  install.packages("devtools")
```

3. Load the devtools package by running the command:

```{r}
  library(devtools)
```

4. Use the function install_github() to download and install the package from the GitHub repository. The syntax for this function is: 

```{r}
  install_github("yolshanskiy/OpenOptions"). 
```

5. Wait for the package to be installed and check the package is installed by running the command installed.packages().

6. Once it is installed, you can load the package by running the command

```{r}
  library(OpenOptions).
```

7. You can now use the functions and data within the package for your analysis.

## Installation via terminal

1. Open the terminal on your computer. Try to skip step 2 and do it only if get error at step 5.

2. Install the R package "devtools" by running the command: 

```
  R -e "install.packages('devtools')"
```

3. Use the command "git clone" to clone the GitHub repository to your local machine. The syntax for this command is:

```
  git clone https://github.com/yolshanskiy/OpenOptions.git.
```  
  
Alternatively, download zip of repository by pressing green code button in the top-right corner.

4. Navigate to the package directory by running the command:

```
  cd path/to/mypackage
```

5. Install the package by running the command

```
  R CMD INSTALL --build OpenOptions.
```
6. Once the package is installed, you can load it in R by running the command

```
  R -e "library(OpenOptions)"
```
