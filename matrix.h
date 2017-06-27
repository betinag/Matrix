#include <cstddef>
#include <vector>

const double MATRIX_DEFAULT_VALUE = 0.0;

class Matrix
{
public:
    typedef std::vector<double>::size_type size_type;

    Matrix();
    Matrix(size_type, size_type, double = MATRIX_DEFAULT_VALUE);
    Matrix(const std::vector<double>&);
    Matrix(const std::vector<std::vector<double> >&);

    Matrix(const Matrix&);
    Matrix(Matrix&&);

    Matrix& operator=(const Matrix&);
    Matrix& operator=(Matrix&&);

    size_type rows() const { return m_rows; }
    size_type cols() const { return m_cols; }

    double& operator()(size_type, size_type);
    double operator()(size_type, size_type) const;

    bool isEmpty() const { return (m_rows == 0 || m_cols == 0); }
    bool isSquare() const { return m_rows == m_cols; }
    bool isFinite() const;

    Matrix operator+(double) const;
    Matrix operator-(double) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const;
    Matrix& operator+=(double);
    Matrix& operator-=(double);
    Matrix& operator*=(double);
    Matrix& operator/=(double);

    Matrix operator-() const;
    friend Matrix operator+(const Matrix&, const Matrix&);
    friend Matrix operator-(const Matrix&, const Matrix&);
    friend Matrix operator*(const Matrix&, const Matrix&);
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);

    Matrix transpose() const;

private:
    size_type m_rows;
    size_type m_cols;
    std::vector<double> m_data;

    int maxRowSize(const std::vector<std::vector<double> >&);
};

inline Matrix operator+(double scalar, const Matrix& rhs)
{
    return rhs + scalar;
}

inline Matrix operator-(double scalar, const Matrix& rhs)
{
    return -rhs + scalar;
}

inline Matrix operator*(double scalar, const Matrix& rhs)
{
    return rhs * scalar;
}
