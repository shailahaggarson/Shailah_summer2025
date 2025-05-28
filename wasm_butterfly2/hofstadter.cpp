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
    
    // For larger matrices, we'll use a basic iterative method
    // This is a simplified implementation - in practice you'd want a more robust solver
    vector<double> eigenvals;
    
    // For tridiagonal matrices, we can use a more direct approach
    // This is a basic implementation that works for the symmetric tridiagonal case
    vector<vector<double>> A = matrix;
    
    // Simple power iteration for dominant eigenvalue, then deflation
    // For the butterfly, we need all eigenvalues, so we'll use a basic QR-like approach
    
    for (int iter = 0; iter < 100; iter++) {
        // Very basic QR step simulation
        for (int i = 0; i < n-1; i++) {
            for (int j = i+1; j < n; j++) {
                if (abs(A[j][i]) > 1e-10) {
                    double c = A[i][i] / sqrt(A[i][i]*A[i][i] + A[j][i]*A[j][i]);
                    double s = A[j][i] / sqrt(A[i][i]*A[i][i] + A[j][i]*A[j][i]);
                    
                    // Apply Givens rotation
                    for (int k = 0; k < n; k++) {
                        double temp1 = c * A[i][k] + s * A[j][k];
                        double temp2 = -s * A[i][k] + c * A[j][k];
                        A[i][k] = temp1;
                        A[j][k] = temp2;
                    }
                    
                    for (int k = 0; k < n; k++) {
                        double temp1 = c * A[k][i] + s * A[k][j];
                        double temp2 = -s * A[k][i] + c * A[k][j];
                        A[k][i] = temp1;
                        A[k][j] = temp2;
                    }
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
        eigs1 = {4.0};
        return {eigs1, eigs2};
    } else if (q == 2) {
        vector<vector<double>> H1 = {{-2, 2}, {2, 2}};
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
        eigs1 = {-4.0};
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
                        Xanti.resize(q);
                        for (int i = 0; i < q; i++) {
                            Xanti[i] = -Xper[q - 1 - i];
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