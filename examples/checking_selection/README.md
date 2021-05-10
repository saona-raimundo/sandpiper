# Motivation
Contrast results with Fyodor's code.
# To do
1. Round

   - Fix parameters
     - population size = 5_000
     - mutation rate = 1.2 e-6
     - beta = 0 (or fixed dominance 0.5)

   - Compute expected polymorphisms ([Heterozygosity](file:///C:/Users/rsaonaur/projects/sandpiper/target/doc/sandpiper/struct.Heterozygosity.html)) for varying values of selection s, with all other parameters fixed. 
      - -10^-2, -10^-3, -10^-4, -10^-5, -10^-6, 0, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2
2. Round

   - Fix parameters
     - population size = 500 
     - mutation rate = 1.2x10-5
     - beta = 0 (or fixed dominance 0.5)
   - Compute expected polymorphisms ([Heterozygosity](file:///C:/Users/rsaonaur/projects/sandpiper/target/doc/sandpiper/struct.Heterozygosity.html)) for varying values of selection s, with all other parameters fixed. 
     - -10^-2, -10^-3, -10^-4, -10^-5, -10^-6, 0, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2
3. Round

   - Fix parameters
     - mu = -0.0000774461, sigma =  0.0000126414, alpha =  0.0
     - population size = 500
     - mutation rate = 1.2x10-5
     - beta = 0 (or fixed dominance 0.5)
   - Compute expected polymorphisms ([Heterozygosity](file:///C:/Users/rsaonaur/projects/sandpiper/target/doc/sandpiper/struct.Heterozygosity.html))
4. Round

   - Fix parameters
     - mu = -0.13458588, sigma =  0.0377358, alpha =  0.0
     - population size = 5000
     - mutation rate = 1.2x10-6
     - beta = 3000
   - Compute expected polymorphisms ([Heterozygosity](file:///C:/Users/rsaonaur/projects/sandpiper/target/doc/sandpiper/struct.Heterozygosity.html))
5. Round

   - Fix parameters
     - population size = 500 
     - mutation rate = 1.2x10-5
     - beta = 3000
   - Compute expected polymorphisms ([Heterozygosity](file:///C:/Users/rsaonaur/projects/sandpiper/target/doc/sandpiper/struct.Heterozygosity.html)) for varying values of selection s, with all other parameters fixed. 
     - -10^-1, -10^-2, -10^-3, -10^-4, -10^-5, -10^-6, -10^-7, 0, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1

# Results

- Numerical results for
  - Mean polymorphisms
- Plots with x-axis in log scale of 
  - Mean