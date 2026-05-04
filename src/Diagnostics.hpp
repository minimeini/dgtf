#pragma once
#include <RcppArmadillo.h>

namespace diagnostics {

inline double effective_sample_size(const arma::vec& draws_in)
{
    arma::vec vals = draws_in.elem(arma::find_finite(draws_in));
    const arma::uword n = vals.n_elem;
    if (n < 2) return static_cast<double>(n);

    vals -= arma::mean(vals);

    arma::uword nfft = 1;
    while (nfft < 2 * n) nfft <<= 1;

    arma::vec padded(nfft, arma::fill::zeros);
    padded.head(n) = vals;

    arma::cx_vec freq = arma::fft(padded);
    freq %= arma::conj(freq);                  // |X(ω)|²
    arma::cx_vec acov_c = arma::ifft(freq);

    arma::vec acov(n);
    for (arma::uword k = 0; k < n; ++k)
        acov(k) = acov_c(k).real() / static_cast<double>(n - k);  // unbiased

    const double gamma0 = acov(0);
    if (!std::isfinite(gamma0) || gamma0 <= 0.0)
        return static_cast<double>(n);

    double sum_rho = 0.0;
    for (arma::uword lag = 1; lag + 1 < n; lag += 2) {
        double pair = acov(lag) / gamma0 + acov(lag + 1) / gamma0;
        if (!std::isfinite(pair) || pair < 0.0) break;             // Geyer IPS
        sum_rho += pair;
    }

    double ess = static_cast<double>(n) / (1.0 + 2.0 * sum_rho);
    return std::min(static_cast<double>(n), std::max(1.0, ess));
}

} // namespace diagnostics
