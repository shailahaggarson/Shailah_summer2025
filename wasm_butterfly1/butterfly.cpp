#include <emscripten/bind.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// Explicitly specify which namespaces we're using
using std::vector;
using std::cos;
using std::sort;
using std::fmod;
using std::sqrt;
using std::pow;
//const double M_PI = 3.14159265358979323846;

// Helper function to compute GCD
int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

// Tridiagonal matrix creation helper
vector<vector<double>> createTridiagonal(int n, const vector<double>& lower, 
                                        const vector<double>& main, 
                                        const vector<double>& upper) {
    vector<vector<double>> matrix(n, vector<double>(n, 0));
    
    for (int i = 0; i < n; i++) {
        if (i > 0) matrix[i][i-1] = lower[i-1];
        matrix[i][i] = main[i];
        if (i < n-1) matrix[i][i+1] = upper[i];
    }
    
    return matrix;
}

// Simplified eigenvalue computation
vector<double> computeEigenvalues(const vector<vector<double>>& matrix) {
    vector<double> eigenvalues;
    if (matrix.empty()) return eigenvalues;
    
    if (matrix.size() == 1) {
        eigenvalues.push_back(matrix[0][0]);
        return eigenvalues;
    }
    
    if (matrix.size() == 2) {
        double a = matrix[0][0], b = matrix[0][1];
        double c = matrix[1][0], d = matrix[1][1];
        double trace = a + d;
        double det = a*d - b*c;
        double disc = trace*trace - 4*det;
        eigenvalues.push_back((trace + sqrt(disc))/2);
        eigenvalues.push_back((trace - sqrt(disc))/2);
        return eigenvalues;
    }
    
    // For larger matrices (simplified)
    for (size_t i = 0; i < matrix.size(); i++) {
        eigenvalues.push_back(matrix[i][i]);
    }
    sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

// Periodic case matrix construction
vector<vector<vector<double>>> Hper(int p, int q, double lambda) {
    if (q == 1) {
        return {{{lambda * 4}}, {}};
    } else if (q == 2) {
        return {{{-2*lambda, 2*lambda}, {2*lambda, 2*lambda}}, {}};
    } else if (q % 2 == 0) {
        int q2 = q / 2;
        vector<double> np(q2+1);
        vector<double> D(q2+1);
        for (int i = 0; i <= q2; i++) {
            np[i] = fmod(i * p, q);
            D[i] = lambda * cos(2 * M_PI * np[i] / q);
        }
        
        vector<double> C = {sqrt(2)};
        for (int i = 0; i < q2-2; i++) C.push_back(1);
        C.push_back(sqrt(2));
        
        auto H1 = createTridiagonal(q2+1, C, D, C);
        vector<vector<double>> H2;
        if (q2 > 1) {
            H2 = vector<vector<double>>(q2-1, vector<double>(q2-1));
            for (int i = 1; i < q2; i++) {
                for (int j = 1; j < q2; j++) {
                    H2[i-1][j-1] = H1[i][j];
                }
            }
        }
        
        return {H1, H2};
    } else {
        int q2 = (q + 1) / 2;
        vector<double> np(q2);
        vector<double> D(q2);
        for (int i = 0; i < q2; i++) {
            np[i] = fmod(i * p, q);
            D[i] = lambda * cos(2 * M_PI * np[i] / q);
        }
        
        vector<double> C = {sqrt(2)};
        for (int i = 0; i < q2-2; i++) C.push_back(1);
        
        auto H1 = createTridiagonal(q2, C, D, C);
        vector<vector<double>> H2;
        if (q2 > 1) {
            H2 = vector<vector<double>>(q2-1, vector<double>(q2-1));
            for (int i = 1; i < q2; i++) {
                for (int j = 1; j < q2; j++) {
                    H2[i-1][j-1] = H1[i][j];
                }
            }
        }
        
        H1[q2-1][q2-1] += 1;
        if (!H2.empty()) {
            H2[q2-2][q2-2] -= 1;
        }
        
        return {H1, H2};
    }
}

// Anti-periodic case matrix construction
vector<vector<vector<double>>> Hanti(int p, int q, double lambda) {
    if (q == 1) {
        return {{{-lambda * 4}}, {}};
    } else if (q == 2) {
        return {{{0}}, {{0}}};
    } else if (q % 2 == 0) {
        int q2 = q / 2;
        vector<double> np(q2);
        vector<double> D(q2);
        for (int i = 0; i < q2; i++) {
            np[i] = fmod((2*i + 1) * p, 2 * q);
            D[i] = lambda * cos(M_PI * np[i] / q);
        }
        
        vector<double> C(q2, 1);
        
        auto H1 = createTridiagonal(q2, C, D, C);
        auto H2 = H1;
        
        H1[0][0] += 1;
        H1[q2-1][q2-1] -= 1;
        H2[0][0] -= 1;
        H2[q2-1][q2-1] += 1;
        
        return {H1, H2};
    } else {
        int q2 = (q + 1) / 2;
        vector<double> np(q2);
        vector<double> D(q2);
        for (int i = 0; i < q2; i++) {
            np[i] = fmod(i * p, q);
            D[i] = -lambda * cos(2 * M_PI * np[i] / q);
        }
        
        vector<double> C = {-sqrt(2)};
        for (int i = 0; i < q2-2; i++) C.push_back(1);
        
        auto H1 = createTridiagonal(q2, C, D, C);
        vector<vector<double>> H2;
        if (q2 > 1) {
            H2 = vector<vector<double>>(q2-1, vector<double>(q2-1));
            for (int i = 1; i < q2; i++) {
                for (int j = 1; j < q2; j++) {
                    H2[i-1][j-1] = H1[i][j];
                }
            }
        }
        
        H1[q2-1][q2-1] -= 1;
        if (!H2.empty()) {
            H2[q2-2][q2-2] += 1;
        }
        
        return {H1, H2};
    }
}

// Calculate spectrum for a given p/q
vector<vector<double>> calculateSpectrum(int p, int q, double lambda) {
    auto result_per = Hper(p, q, lambda);
    auto H1per = result_per[0];
    auto H2per = result_per[1];
    
    vector<double> Xper;
    
    auto eig1 = computeEigenvalues(H1per);
    Xper.insert(Xper.end(), eig1.begin(), eig1.end());
    
    if (!H2per.empty()) {
        auto eig2 = computeEigenvalues(H2per);
        Xper.insert(Xper.end(), eig2.begin(), eig2.end());
    }
    
    sort(Xper.begin(), Xper.end());
    
    vector<double> Xanti;
    if (q % 2 == 0) {
        auto result_anti = Hanti(p, q, lambda);
        auto H1anti = result_anti[0];
        auto H2anti = result_anti[1];
        
        auto eig1 = computeEigenvalues(H1anti);
        Xanti.insert(Xanti.end(), eig1.begin(), eig1.end());
        
        if (!H2anti.empty()) {
            auto eig2 = computeEigenvalues(H2anti);
            Xanti.insert(Xanti.end(), eig2.begin(), eig2.end());
        }
        
        sort(Xanti.begin(), Xanti.end());
    } else {
        Xanti.resize(Xper.size());
        for (size_t i = 0; i < Xper.size(); i++) {
            Xanti[i] = -Xper[Xper.size() - 1 - i];
        }
    }
    
    return {Xper, Xanti};
}

// Explicitly qualify the EMSCRIPTEN_BINDINGS block
EMSCRIPTEN_BINDINGS(butterfly_module) {
    emscripten::function("calculateSpectrum", &calculateSpectrum);
    emscripten::function("gcd", &gcd);
    
    emscripten::register_vector<double>("VectorDouble");
    emscripten::register_vector<vector<double>>("VectorVectorDouble");
}