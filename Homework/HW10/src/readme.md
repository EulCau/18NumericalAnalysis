# Experiment 10

## Adams-Bashforth Method (AB5)

To run the experiment using the 5th-order Adams-Bashforth method, run the following file: `ab5/ab5_test.m`

### Initial Condition:
- ODE: $ x' = \frac{t - e^{-t}}{x + e^x} $
- Initial value: `t0 = 0`, `x0 = 1`
- Reference solution at `t = 1` is computed via the implicit relation:
  $ x^2 - t^2 + 2e^x - 2e^{-t} + 1 - 2e $

## Adaptive RKF Method

To run the adaptive Runge-Kutta-Fehlberg method (RKF), run the following file: `rkf/rkf_text.m`

### Initial Condition:
- ODE: $ y' = e^{yx} - \cos(y - x) $
- Initial value: `x0 = 1`, `y0 = 3`
- The solver supports interactive interpolation after computation.
