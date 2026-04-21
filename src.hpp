#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <vector>
#include <utility>

class matrix {
private:
    int m, n;
    fraction **data;

    void allocate(int rows, int cols) {
        m = rows; n = cols;
        if (m <= 0 || n <= 0) { data = nullptr; return; }
        data = new fraction*[m];
        try {
            for (int i = 0; i < m; ++i) data[i] = new fraction[n];
        } catch (...) {
            if (data) {
                for (int i = 0; i < m; ++i) if (data[i]) delete [] data[i];
                delete [] data;
            }
            data = nullptr; m = n = 0;
            throw;
        }
    }

public:
    matrix() { m = n = 0; data = nullptr; }

    // 构造 m_*n_ 的矩阵，元素设为0
    matrix(int m_, int n_) { allocate(m_, n_); }

    // 拷贝构造
    matrix(const matrix &obj) {
        allocate(obj.m, obj.n);
        if (data) {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    data[i][j] = obj.data[i][j];
        }
    }

    // 移动拷贝构造
    matrix(matrix &&obj) noexcept { m = obj.m; n = obj.n; data = obj.data; obj.m = obj.n = 0; obj.data = nullptr; }

    // 析构
    ~matrix() {
        if (data) {
            for (int i = 0; i < m; ++i) delete [] data[i];
            delete [] data;
        }
        data = nullptr; m = n = 0;
    }

    // 赋值
    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        matrix tmp(obj);
        std::swap(m, tmp.m); std::swap(n, tmp.n); std::swap(data, tmp.data);
        return *this;
    }

    // 访问(i,j): i 1-based, j 0-based
    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) throw matrix_error();
        return data[i-1][j];
    }

    const fraction &at(int i, int j) const {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) throw matrix_error();
        return data[i-1][j];
    }

    int rows() const { return m; }
    int cols() const { return n; }

    // 乘法
    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n == 0 || rhs.m == 0) throw matrix_error();
        if (lhs.n != rhs.m) throw matrix_error();
        matrix res(lhs.m, rhs.n);
        for (int i = 1; i <= lhs.m; ++i) {
            for (int k = 0; k < lhs.n; ++k) {
                fraction a = lhs.at(i, k);
                if (a == fraction(0)) continue;
                for (int j = 0; j < rhs.n; ++j) {
                    res(i, j) = res(i, j) + a * rhs.at(k+1, j);
                }
            }
        }
        return res;
    }

    // 转置
    matrix transposition() {
        if (m == 0 || n == 0 || data == nullptr) throw matrix_error();
        matrix t(n, m);
        for (int i = 1; i <= m; ++i) {
            for (int j = 0; j < n; ++j) {
                t(j+1, i-1) = at(i, j);
            }
        }
        return t;
    }

    // 行列式 (高斯消元)
    fraction determination() {
        if (m == 0 || n == 0 || data == nullptr) throw matrix_error();
        if (m != n) throw matrix_error();
        matrix a(*this);
        fraction det(1);
        int sign = 1;
        for (int col = 0, row = 0; col < n && row < m; ++col, ++row) {
            int piv = -1;
            for (int i = row; i < m; ++i) {
                if (!(a.data[i][col] == fraction(0))) { piv = i; break; }
            }
            if (piv == -1) return fraction(0);
            if (piv != row) { std::swap(a.data[piv], a.data[row]); sign = -sign; }
            fraction pivot = a.data[row][col];
            // eliminate below
            for (int i = row+1; i < m; ++i) {
                if (a.data[i][col] == fraction(0)) continue;
                fraction factor = a.data[i][col] / pivot;
                for (int j = col; j < n; ++j) a.data[i][j] = a.data[i][j] - factor * a.data[row][j];
            }
            det = det * pivot;
        }
        if (sign < 0) det = fraction(0) - det;
        return det;
    }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction; // A, C

    void build_matrices(int from[], int to[], fraction resistance[]) {
        adjacency = matrix(connection_size, interface_size); // m x n
        conduction = matrix(connection_size, connection_size); // m x m diag
        for (int e = 0; e < connection_size; ++e) {
            int u = from[e] - 1;
            int v = to[e] - 1;
            adjacency(e+1, u) = fraction(1);
            adjacency(e+1, v) = fraction(0) - fraction(1);
            conduction(e+1, e) = fraction(1) / resistance[e];
        }
    }

    matrix laplacian() {
        matrix AT = adjacency.transposition();
        matrix CA = conduction * adjacency;
        matrix G = AT * CA; // n x n
        return G;
    }

    static std::vector<fraction> solve_linear(matrix M, std::vector<fraction> b) {
        int n = M.rows();
        if (n != M.cols() || (int)b.size() != n) throw matrix_error();
        // Gauss-Jordan
        int row = 0;
        for (int col = 0; col < n && row < n; ++col) {
            int piv = -1;
            for (int i = row; i < n; ++i) if (!(M(i+1, col) == fraction(0))) { piv = i; break; }
            if (piv == -1) continue; // dependent column
            if (piv != row) {
                for (int j = 0; j < n; ++j) { fraction tmp = M(row+1, j); M(row+1, j) = M(piv+1, j); M(piv+1, j) = tmp; }
                std::swap(b[row], b[piv]);
            }
            fraction pivot = M(row+1, col);
            for (int j = col; j < n; ++j) M(row+1, j) = M(row+1, j) / pivot;
            b[row] = b[row] / pivot;
            for (int i = 0; i < n; ++i) if (i != row) {
                if (M(i+1, col) == fraction(0)) continue;
                fraction factor = M(i+1, col);
                for (int j = col; j < n; ++j) M(i+1, j) = M(i+1, j) - factor * M(row+1, j);
                b[i] = b[i] - factor * b[row];
            }
            ++row;
        }
        return b;
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[])
        : interface_size(interface_size_), connection_size(connection_size_), adjacency(), conduction() {
        build_matrices(from, to, resistance);
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1 - 1;
        int b = interface_id2 - 1;
        if (a == b) return fraction(0);
        matrix G = laplacian();
        int n = interface_size;
        int k = n - 1; // ground last node
        matrix R(k, k);
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                R(i+1, j) = G.at(i+1, j);
        std::vector<fraction> I(k, fraction(0));
        if (a < n-1) I[a] = I[a] + fraction(1);
        if (b < n-1) I[b] = I[b] - fraction(1);
        std::vector<fraction> U = solve_linear(R, I);
        fraction ua = (a < n-1) ? U[a] : fraction(0);
        fraction ub = (b < n-1) ? U[b] : fraction(0);
        return ua - ub;
    }

    fraction get_voltage(int id, fraction current[]) {
        matrix G = laplacian();
        int n = interface_size;
        int k = n - 1;
        matrix R(k, k);
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                R(i+1, j) = G.at(i+1, j);
        std::vector<fraction> I(k);
        for (int i = 0; i < k; ++i) I[i] = current[i]; // ground last node current omitted
        std::vector<fraction> U = solve_linear(R, I);
        return U[id-1]; // id < n ensured
    }

    fraction get_power(fraction voltage[]) {
        // P = U^T G U
        matrix G = laplacian();
        int n = interface_size;
        std::vector<fraction> y(n, fraction(0));
        for (int i = 1; i <= n; ++i) {
            fraction s(0);
            for (int j = 0; j < n; ++j) s = s + G.at(i, j) * voltage[j];
            y[i-1] = s;
        }
        fraction total(0);
        for (int i = 0; i < n; ++i) total = total + voltage[i] * y[i];
        return total;
    }
};

#endif // SRC_HPP
