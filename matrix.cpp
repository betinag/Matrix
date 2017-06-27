#include <cmath>
#include <stdexcept>
#include "matrix.h"

Matrix::Matrix()
    : m_rows(0)
    , m_cols(0)
{
}

Matrix::Matrix(size_type rows, size_type cols, double value)
    : m_rows(rows)
    , m_cols(cols)
    , m_data(rows * cols, value)
{
}

Matrix::Matrix(const std::vector<double>& data)
    : m_rows(0)
    , m_cols(0)
{
    if (!data.empty())
    {
        m_rows = 1;
        m_cols = data.size();
        m_data = data;
    }
}

Matrix::Matrix(const std::vector<std::vector<double> >& data)
    : m_rows(0)
    , m_cols(0)
{
    if (!data.empty())
    {
        m_rows = data.size();
        m_cols = maxRowSize(data);
        m_data.reserve(m_rows * m_cols);
        for (size_type i = 0; i < m_rows; ++i)
        {
            m_data.insert(m_data.end(), data[i].begin(), data[i].end());
            m_data.insert(m_data.end(), m_cols - data[i].size(), MATRIX_DEFAULT_VALUE);
        }
    }
}

int Matrix::maxRowSize(const std::vector<std::vector<double> >& data)
{
    size_type maxSize = 0;
    size_type rows = data.size();
    for (size_type i = 0; i < rows; ++i)
    {
        maxSize = std::max(data[i].size(), maxSize);
    }

    return maxSize;
}

Matrix::Matrix(const Matrix& rhs)
    : m_rows(rhs.m_rows)
    , m_cols(rhs.m_cols)
    , m_data(rhs.m_data)
{
}

Matrix::Matrix(Matrix&& rhs)
    : m_rows(rhs.m_rows)
    , m_cols(rhs.m_cols)
    , m_data(std::move(rhs.m_data))
{
    rhs.m_rows = 0;
    rhs.m_cols = 0;
}

Matrix& Matrix::operator=(const Matrix& rhs)
{
    if (this != &rhs)
    {
        m_rows = rhs.m_rows;
        m_cols = rhs.m_cols;
        m_data = rhs.m_data;
    }

    return *this;
}

Matrix& Matrix::operator=(Matrix&& rhs)
{
    m_rows = rhs.m_rows;
    m_cols = rhs.m_cols;
    m_data = std::move(rhs.m_data);

    rhs.m_rows = 0;
    rhs.m_cols = 0;

    return *this;
}

double Matrix::operator()(size_type row, size_type col) const
{
    if (row >= m_rows || col >= m_cols)
    {
        throw std::out_of_range("Trying to access an element outside the bounds of the matrix.");
    }

    return m_data[row * m_cols + col];
}

double& Matrix::operator()(size_type row, size_type col)
{
    if (row >= m_rows || col >= m_cols)
    {
        throw std::out_of_range("Trying to access an element outside the bounds of the matrix.");
    }

    return m_data[row * m_cols + col];
}

bool Matrix::isFinite() const
{
    size_type elements = m_data.size();
    for(size_type i = 0; i < elements; ++i)
    {
        if(!std::isfinite(m_data[i]))
        {
            return false;
        }
    }

    return true;
}

Matrix Matrix::operator+(double scalar) const
{
    Matrix result = *this;
    result += scalar;

    return result;
}

Matrix Matrix::operator-(double scalar) const
{
    return (*this) + (-1 * scalar);
}

Matrix Matrix::operator*(double scalar) const
{
    Matrix result = *this;
    result *= scalar;

    return result;
}

Matrix Matrix::operator/(double scalar) const
{
    Matrix result = *this;
    result /= scalar;

    return result;
}

Matrix& Matrix::operator+=(double scalar)
{
    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] += scalar;
    }

    return *this;
}

Matrix& Matrix::operator-=(double scalar)
{
    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] -= scalar;
    }

    return *this;
}

Matrix& Matrix::operator*=(double scalar)
{
    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] *= scalar;
    }

    return *this;
}

Matrix& Matrix::operator/=(double scalar)
{
    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] /= scalar;
    }

    return *this;
}

Matrix Matrix::operator-() const
{
    return (*this) * (-1);
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
    if (lhs.m_rows != rhs.m_rows || lhs.m_cols != rhs.m_cols)
    {
        throw std::invalid_argument("Incompatible matrix dimensions.");
    }

    Matrix result(lhs.m_rows, lhs.m_cols);
    Matrix::size_type elements = lhs.m_data.size();
    for (Matrix::size_type i = 0; i < elements; ++i)
    {
        result.m_data[i] = lhs.m_data[i] + rhs.m_data[i];
    }

    return result;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
    if (lhs.m_rows != rhs.m_rows || lhs.m_cols != rhs.m_cols)
    {
        throw std::invalid_argument("Incompatible matrix dimensions.");
    }

    Matrix result(lhs.m_rows, lhs.m_cols);
    Matrix::size_type elements = lhs.m_data.size();
    for (Matrix::size_type i = 0; i < elements; ++i)
    {
        result.m_data[i] = lhs.m_data[i] - rhs.m_data[i];
    }

    return result;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
    if (lhs.m_cols != rhs.m_rows)
    {
        throw std::invalid_argument("Incompatible matrix dimensions.");
    }

    Matrix result(lhs.m_rows, rhs.m_cols);
    for (Matrix::size_type i = 0; i < lhs.m_rows; ++i)
    {
        for (Matrix::size_type j = 0; j < rhs.m_cols; ++j)
        {
            result(i, j) = 0.0;
            for (Matrix::size_type k = 0; k < lhs.m_cols; ++k)
                result(i, j) += lhs(i, k) * rhs(k, j);
        }
    }

    return result;
}

Matrix& Matrix::operator+=(const Matrix& rhs)
{
    if (m_rows != rhs.m_rows || m_cols != rhs.m_cols)
    {
        throw std::invalid_argument("Incompatible matrix dimensions.");
    }

    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] += rhs.m_data[i];
    }

    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs)
{
    if (m_rows != rhs.m_rows || m_cols != rhs.m_cols)
    {
        throw std::invalid_argument("Incompatible matrix dimensions.");
    }

    size_type elements = m_data.size();
    for (size_type i = 0; i < elements; ++i)
    {
        m_data[i] -= rhs.m_data[i];
    }

    return *this;
}

Matrix Matrix::transpose() const
{
    if (isEmpty())
    {
        return *this;
    }

    Matrix transposed(m_cols, m_rows);
    for (size_type i = 0; i < m_rows; ++i)
    {
        for (size_type j = 0; j < m_cols; ++j)
        {
            transposed(j, i) = this->operator()(i, j);
        }
    }

    return transposed;
}
