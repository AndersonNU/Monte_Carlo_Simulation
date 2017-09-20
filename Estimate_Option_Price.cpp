// Using antithetic variance method to estimate a European option price
// Running time is 58 seconds
#include "Declarations.h"

main() {

   double r, t, T, mu, sigma, dt, sqrtdt, S, Sa, S0, V, Vbar, V2bar,
          elapsed_time, t_star, stdhatV, error, epsilon, n, Discount_factor,
          U, W, Wa;
   double C,Ca, K;
   int i, N, done, test;

   K = 110.0;
   // Time to expiration.
   T = 0.5;

   // Number of stock price periods.
   N = 50;

   // Time increment per period.
   dt = T / N;

   // Compute this oft-used value once and for all.
   sqrtdt = sqrt(dt);

   // Risk-free interest rate.
   r = 0.05;

   // Compute the oft-used discount factor once and for all.
   Discount_factor = exp(-r*T);

   // Stock price volatility.
   sigma = .30;

   // Drift term.
   mu = r - log(cosh(sigma*sqrtdt)) / dt;

   // Initial stock price.
   S0 = 100.0;

   // Specify the 95% error tolerance.
   epsilon = 0.005;

   // Start the clock to time the computations.
   Time ();

   // Seed the RNG.
   MTUniform (1);

   // Print column headings for output to execution window.
   printf ("         n       Vbar        +/-        t       t*\n");

   // Initialize certain values.
   V2bar = Vbar = done = n = test = 0;

   // Begin the simulation loop.
   while (!done) {

      // Initialize the stock price.
      S = Sa = S0;

      // Initialize the random walk.
      W = Wa = 0;

      // Initialize time.
      t = 0;

      // Simulate a stock price path.  Go forward period-by-period computing
      //   the next stock price.
      for (i = 1; i <= N; i++) {

         // Advance the random walk.
         U = MTUniform (0);
         W += (U <= 0.5 ? -1 : +1);
         Wa = -W ;
         // Advance time by one period.
         t += dt;
                     // Compute the next stock price.
         S = S0 * exp (mu*t + sigma*sqrtdt*W);
         Sa =S0 * exp (mu*t + sigma*sqrtdt*Wa);

      }
      // S now is the value of the stock at time T for the simulated price path.




         // Compute the payout of the option
         C = (S>= K ? (S-K) : 0);
         Ca= (Sa>= K ? (Sa-K) : 0);
      // Discount back to time 0.
      V  = Discount_factor * ((C+Ca)/2);

      // Update the simulation counter and the sample moments of V at time T.
      n ++;
      Vbar  = ((n-1)*Vbar + V) / n;
      V2bar = ((n-1)*V2bar + V*V) / n;

      // Update the error condition test counter.
      test ++;

      // Test the error condition periodically.
      if (test == 100000) {

         // Estimate the standard deviation of S.
         stdhatV  = sqrt (V2bar - Vbar*Vbar);

          // Estimate the error of the Vbar estimator.
         error = 1.96 * stdhatV / sqrt(n);

         // Compute the elapsed time and estimate the time of completion.
         elapsed_time = Time ();
         t_star = elapsed_time * pow (error / epsilon, 2);

         // Report.
         printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);

         // Reset the "test" counter and see if error tolerance is met.
         test = 0;
         if (error < epsilon) {
            done = 1;
         }

      }

   }

   Exit ();

}

#include "Definitions.h"

