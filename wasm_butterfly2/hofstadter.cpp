#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>

using namespace std;

// Structure to hold matrix data
struct SparseMatrix {
    vector<double> data;
    vector<int> row_indices;
    vector<int> col_indices;
    int rows, cols;
    
    SparseMatrix(int r, int c) : rows(r), cols(c) {}
    
    void addElement(int row, int col, double value) {
        if (abs(value) > 1e-12) {  // Only store non-zero elements
            data.push_back(value);
            row_indices.push_back(row);
            col_indices.push_back(col);
        }
    }
    
    // Convert to dense matrix for eigenvalue computation
    vector<vector<double>> toDense() const {
        vector<vector<double>> dense(rows, vector<double>(cols, 0.0));
        for (size_t i = 0; i < data.size(); i++) {
            dense[row_indices[i]][col_indices[i]] = data[i];
        }
        return dense;
    }
};

// GCD function
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Simple eigenvalue computation using QR algorithm (simplified version)
vector<double> computeEigenvalues(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    if (n == 0) return {};
    
    // For small matrices, we can use a simplified approach
    if (n == 1) {
        return {matrix[0][0]};
    }
    if (n == 2) {
        // For 2x2 matrix, compute eigenvalues directly
        double a = matrix[0][0], b = matrix[0][1];
        double c = matrix[1][0], d = matrix[1][1];
        double trace = a + d;
        double det = a*d - b*c;
        double discriminant = trace*trace - 4*det;
        if (discriminant < 0) {
            // Shouldn't happen for symmetric matrices
            return {trace/2, trace/2};
        }
        double sqrt_disc = sqrt(discriminant);
        return {(trace + sqrt_disc)/2, (trace - sqrt_disc)/2};
    }
    
    // For larger matrices, we'll use a basic iterative method
    vector<double> eigenvals;
    vector<vector<double>> A = matrix;
    
    // Simple QR iteration (very basic implementation)
    for (int iter = 0; iter < 100; iter++) {
        // QR decomposition (simplified)
        vector<vector<double>> Q(n, vector<double>(n, 0.0));
        vector<vector<double>> R(n, vector<double>(n, 0.0));
        
        // Gram-Schmidt process (simplified)
        for (int j = 0; j < n; j++) {
            vector<double> v = A[j];
            for (int k = 0; k < j; k++) {
                double dot = 0.0;
                for (int i = 0; i < n; i++) {
                    dot += Q[i][k] * A[j][i];
                }
                for (int i = 0; i < n; i++) {
                    v[i] -= dot * Q[i][k];
                }
            }
            double norm = 0.0;
            for (int i = 0; i < n; i++) {
                norm += v[i] * v[i];
            }
            norm = sqrt(norm);
            if (norm > 1e-12) {
                for (int i = 0; i < n; i++) {
                    Q[i][j] = v[i] / norm;
                }
            }
        }
        
        // Compute R = Q^T * A
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                R[i][j] = 0.0;
                for (int k = 0; k < n; k++) {
                    R[i][j] += Q[k][i] * A[k][j];
                }
            }
        }
        
        // Update A = R * Q
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = 0.0;
                for (int k = 0; k < n; k++) {
                    A[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
    }
    
    // Extract diagonal elements as eigenvalue approximations
    for (int i = 0; i < n; i++) {
        eigenvals.push_back(A[i][i]);
    }
    
    sort(eigenvals.begin(), eigenvals.end());
    return eigenvals;
}

// HPER function - Periodic H matrix
pair<vector<double>, vector<double>> Hper(int p, int q, double lambda) {
    vector<double> eigs1, eigs2;
    
    if (q == 1) {
        eigs1 = {2.0 + lambda};  // Fixed from original MATLAB code
        return {eigs1, eigs2};
    } else if (q == 2) {
        vector<vector<double>> H1 = {{-lambda, 2.0}, {2.0, lambda}};
        eigs1 = computeEigenvalues(H1);
        return {eigs1, eigs2};
    } else if (q % 2 == 0) {
        // Even q case
        int q2 = q / 2;
        SparseMatrix H1(q2 + 1, q2 + 1);
        
        // Compute diagonal elements
        for (int i = 0; i <= q2; i++) {
            int np = (i * p) % q;
            double diag_val = lambda * cos(2.0 * M_PI * np / q);
            H1.addElement(i, i, diag_val);
        }
        
        // Off-diagonal elements
        for (int i = 0; i < q2; i++) {
            double off_diag = (i == 0 || i == q2 - 1) ? sqrt(2.0) : 1.0;
            H1.addElement(i, i + 1, off_diag);
            H1.addElement(i + 1, i, off_diag);
        }
        
        auto dense1 = H1.toDense();
        eigs1 = computeEigenvalues(dense1);
        
        // H2 is submatrix
        if (q2 > 1) {
            vector<vector<double>> H2_dense(q2 - 1, vector<double>(q2 - 1));
            for (int i = 0; i < q2 - 1; i++) {
                for (int j = 0; j < q2 - 1; j++) {
                    H2_dense[i][j] = dense1[i + 1][j + 1];
                }
            }
            eigs2 = computeEigenvalues(H2_dense);
        }
    } else {
        // Odd q case
        int q2 = (q + 1) / 2;
        SparseMatrix H1(q2, q2);
        
        // Diagonal elements
        for (int i = 0; i < q2; i++) {
            int np = (i * p) % q;
            double diag_val = lambda * cos(2.0 * M_PI * np / q);
            H1.addElement(i, i, diag_val);
        }
        
        // Off-diagonal elements
        for (int i = 0; i < q2 - 1; i++) {
            double off_diag = (i == 0) ? sqrt(2.0) : 1.0;
            H1.addElement(i, i + 1, off_diag);
            H1.addElement(i + 1, i, off_diag);
        }
        
        auto dense1 = H1.toDense();
        dense1[q2-1][q2-1] += 1.0;  // H1(q2,q2) = H1(q2,q2) + 1
        eigs1 = computeEigenvalues(dense1);
        
        // H2 is submatrix with modification
        if (q2 > 1) {
            vector<vector<double>> H2_dense(q2 - 1, vector<double>(q2 - 1));
            for (int i = 0; i < q2 - 1; i++) {
                for (int j = 0; j < q2 - 1; j++) {
                    H2_dense[i][j] = dense1[i + 1][j + 1];
                }
            }
            H2_dense[q2-2][q2-2] -= 1.0;  // H2(q2-1,q2-1) = H2(q2-1,q2-1) - 1
            eigs2 = computeEigenvalues(H2_dense);
        }
    }
    
    return {eigs1, eigs2};
}

// HANTI function - Anti-periodic H matrix
pair<vector<double>, vector<double>> Hanti(int p, int q, double lambda) {
    vector<double> eigs1, eigs2;
    
    if (q == 1) {
        eigs1 = {-lambda - 2.0};  // Fixed from original MATLAB code
        return {eigs1, eigs2};
    } else if (q == 2) {
        eigs1 = {0.0};
        eigs2 = {0.0};
        return {eigs1, eigs2};
    } else if (q % 2 == 0) {
        // Even q case
        int q2 = q / 2;
        SparseMatrix H1(q2, q2), H2(q2, q2);
        
        // Diagonal elements
        for (int i = 0; i < q2; i++) {
            int np = ((i + 1) * 2 * p) % (2 * q);
            double diag_val = lambda * cos(M_PI * np / q);
            H1.addElement(i, i, diag_val);
            H2.addElement(i, i, diag_val);
        }
        
        // Off-diagonal elements
        for (int i = 0; i < q2 - 1; i++) {
            H1.addElement(i, i + 1, 1.0);
            H1.addElement(i + 1, i, 1.0);
            H2.addElement(i, i + 1, 1.0);
            H2.addElement(i + 1, i, 1.0);
        }
        
        auto dense1 = H1.toDense();
        auto dense2 = H2.toDense();
        
        // Boundary modifications
        dense1[0][0] += 1.0;
        dense1[q2-1][q2-1] -= 1.0;
        dense2[0][0] -= 1.0;
        dense2[q2-1][q2-1] += 1.0;
        
        eigs1 = computeEigenvalues(dense1);
        eigs2 = computeEigenvalues(dense2);
    } else {
        // Odd q case - similar to Hper but with sign changes
        int q2 = (q + 1) / 2;
        SparseMatrix H1(q2, q2);
        
        for (int i = 0; i < q2; i++) {
            int np = (i * p) % q;
            double diag_val = -lambda * cos(2.0 * M_PI * np / q);
            H1.addElement(i, i, diag_val);
        }
        
        for (int i = 0; i < q2 - 1; i++) {
            double off_diag = (i == 0) ? -sqrt(2.0) : 1.0;
            H1.addElement(i, i + 1, off_diag);
            H1.addElement(i + 1, i, off_diag);
        }
        
        auto dense1 = H1.toDense();
        dense1[q2-1][q2-1] -= 1.0;
        eigs1 = computeEigenvalues(dense1);
        
        if (q2 > 1) {
            vector<vector<double>> H2_dense(q2 - 1, vector<double>(q2 - 1));
            for (int i = 0; i < q2 - 1; i++) {
                for (int j = 0; j < q2 - 1; j++) {
                    H2_dense[i][j] = dense1[i + 1][j + 1];
                }
            }
            H2_dense[q2-2][q2-2] += 1.0;
            eigs2 = computeEigenvalues(H2_dense);
        }
    }
    
    return {eigs1, eigs2};
}

// Main computation function that returns all spectrum data
extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double* computeSpectrum(int max_q, double lambda, int* result_size) {
        vector<double> all_data;
        
        // Add zero line
        auto [eigs1_0, eigs2_0] = Hper(0, 1, lambda);
        all_data.push_back(0.0);  // theta = 0
        all_data.push_back(static_cast<double>(eigs1_0.size()));
        for (double eig : eigs1_0) {
            all_data.push_back(eig);
        }
        
        // Process all fractions p/q
        for (int q = 1; q <= max_q; q++) {
            for (int p = 1; p < q; p++) {
                if (gcd(p, q) == 1) {
                    double theta = static_cast<double>(p) / q;
                    
                    auto [Xper1, Xper2] = Hper(p, q, lambda);
                    vector<double> Xper = Xper1;
                    Xper.insert(Xper.end(), Xper2.begin(), Xper2.end());
                    sort(Xper.begin(), Xper.end());
                    
                    vector<double> Xanti;
                    if (q % 2 == 0) {
                        auto [Xanti1, Xanti2] = Hanti(p, q, lambda);
                        Xanti = Xanti1;
                        Xanti.insert(Xanti.end(), Xanti2.begin(), Xanti2.end());
                        sort(Xanti.begin(), Xanti.end());
                    } else {
                        // Xanti = -Xper(q:-1:1) in MATLAB notation
                        Xanti.resize(Xper.size());
                        for (size_t i = 0; i < Xper.size(); i++) {
                            Xanti[i] = -Xper[Xper.size() - 1 - i];
                        }
                    }
                    
                    // Store data: theta, count, eigenvalues
                    all_data.push_back(theta);
                    all_data.push_back(static_cast<double>(Xper.size() + Xanti.size()));
                    for (double x : Xper) all_data.push_back(x);
                    for (double x : Xanti) all_data.push_back(x);
                }
            }
        }
        
        *result_size = all_data.size();
        
        // Allocate and copy data
        double* result = (double*)malloc(all_data.size() * sizeof(double));
        copy(all_data.begin(), all_data.end(), result);
        
        return result;
    }
    
    EMSCRIPTEN_KEEPALIVE
    void freeMemory(double* ptr) {
        free(ptr);
    }
}