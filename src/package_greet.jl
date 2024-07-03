function package_greet()

    greeting = "
                   Hello!
                   Welcome to the SigTools package for Julia. A package optimized for Signal Processing and Data Analysis.

                   version 0.0.1
                   Developed by Sam Kramer
                   March 2024
               "

    println(greeting)
    return nothing
end

function package_list()

    list = "
               This package includes functions for:
                   Convolutions -- Both 1D and 2D
                   Correlations -- Both 1D and 2D
                       - Normalized
                       - Biased
                       - Unbiased
                       - None
                   Kernel Functions
                       - Linear --> Absolute Value, Similarity Score
                       - Nonlinear --> Polynomial, Gaussian

           "

    println(list)
    return nothing
end

