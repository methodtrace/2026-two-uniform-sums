/* Step inequalities I(u_J) on the grid 1/K at a lattice law, exhaustive over all step sets J of {1, ..., K-1}.
   The law has masses w_m / D at the points (t_m + 1) / K, t_m in {0, ..., 2K-2}; for u_J(i) = #{j in J : j <= i}
   and v(t) = min_{i+j=t} u_J(i) + u_J(j) the inequality I(u_J) reads K sum_m w_m v(t_m) <= 2 D sum_i u_J(i).
   Exact 64-bit integer arithmetic. Prints one JSON line with the number of step sets, the number violated and the
   largest value of K sum_m w_m v(t_m) - 2 D sum_i u_J(i) (with the step set attaining it, as a bit mask).
   usage: grid_steps K D t_1 w_1 t_2 w_2 ...        build: cc -O2 -o grid_steps grid_steps.c */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
  if (argc < 5 || (argc - 3) % 2) { fprintf(stderr, "usage: grid_steps K D t1 w1 [t2 w2 ...]\n"); return 2; }
  int K = atoi(argv[1]), D = atoi(argv[2]), n = (argc - 3) / 2, t[32], w[32], u[64];
  if (K < 2 || K > 40 || n > 32) { fprintf(stderr, "K must be in 2..40 and at most 32 atoms\n"); return 2; }
  for (int m = 0; m < n; m++) { t[m] = atoi(argv[3 + 2 * m]); w[m] = atoi(argv[4 + 2 * m]); }
  long long best = -(1LL << 60), cnt = 0, viol = 0; unsigned long long bestJ = 0;
  for (unsigned long long J = 0; J < (1ULL << (K - 1)); J++) {
    int c = 0; long long su = 0;
    for (int i = 0; i < K; i++) { if (i >= 1 && ((J >> (i - 1)) & 1ULL)) c++; u[i] = c; su += c; }
    long long lhs = 0;
    for (int m = 0; m < n; m++) {
      int tt = t[m], lo = tt - (K - 1) > 0 ? tt - (K - 1) : 0, hi = tt < K - 1 ? tt : K - 1, mn = 1 << 30;
      for (int i = lo; i <= hi; i++) { int s = u[i] + u[tt - i]; if (s < mn) mn = s; }
      lhs += (long long)w[m] * mn;
    }
    long long diff = lhs * K - 2LL * D * su;
    cnt++; if (diff > 0) viol++;
    if (diff > best) { best = diff; bestJ = J; }
  }
  printf("{\"K\": %d, \"step_sets\": %lld, \"violated\": %lld, \"max_lhs_minus_rhs\": %lld, \"denominator\": %lld, \"argmax_mask\": %llu}\n",
         K, cnt, viol, best, (long long)D * K, bestJ);
  return 0;
}
