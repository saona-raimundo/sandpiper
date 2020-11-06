# Checking

If s = infinity, then we know that allele frequency must be 1 and therefore the polymorphism is 0.

If s = - infinity, then allele frequency must be 0 and therefore the polymorphism is 0 too.

If s = 0, then allele frequency is the most random and concentrates around 1/2, therefore the polymorphism is positive.

I will check the transition first, by fixing s, and then try to simulate this with some parameters of the model.

# Next steps

-  alpha: 0, -2, -4
  beta: 0, 10, 50, 100, 1000, 5000
  
  with each of these:
  
  - [ ] [Running]
    mu: -0.2000000000, -0.1523191619, -0.1160056354, -0.0883494058, -0.0672865373, -0.0512451448, -0.0390280876, -0.0297236279, -0.0226373905, -0.0172405417, -0.0131303243, -0.0100000000
    sigma: 0.0100000000, 0.0131303243, 0.0172405417, 0.0226373905, 0.0297236279, 0.0390280876, 0.0512451448, 0.0672865373, 0.0883494058, 0.1160056354, 0.1523191619, 0.2000000000
  - [ ] mu: -0.0100000000, -0.0056838040, -0.0032305627, -0.0018361885, -0.0010436536, -0.0005931922, -0.0003371588, -0.0001916345, -0.0001089213, -0.0000619087, -0.0000351877, -0.0000200000 
    sigma: 0.0100000000, 0.0131303243, 0.0172405417, 0.0226373905, 0.0297236279, 0.0390280876, 0.0512451448, 0.0672865373, 0.0883494058, 0.1160056354, 0.1523191619, 0.2000000000
  - [ ] mu: -0.2000000000, -0.1523191619, -0.1160056354, -0.0883494058, -0.0672865373, -0.0512451448, -0.0390280876, -0.0297236279, -0.0226373905, -0.0172405417, -0.0131303243, -0.0100000000
    sigma: 0.0000100000, 0.0000187382, 0.0000351119, 0.0000657933, 0.0001232847, 0.0002310130, 0.0004328761, 0.0008111308, 0.0015199111, 0.0028480359, 0.0053366992, 0.0100000000

# To do

- PR for rand_distr for a SkewNormal
- Implement in R64
- How do I even get a NaN??
- Test correctness of SkewNormal simulation
- Add StandardNormal as a struct?
- Convert Parameters to an Array1 in the UnitSphere