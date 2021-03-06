#ifndef COMPLEX_ARITHMETIC_H
#define COMPLEX_ARITHMETIC_H

#include <cmath>
#include <vector>
#include <sstream>

class ComplexNumber
{
public:
    ComplexNumber(double x, double y) :
            _x(x),
            _y(y)
    {}

    ComplexNumber() :
            _x(std::nan("0")),
            _y(std::nan("0"))
    {}

    double re() const { return _x; }
    double im() const { return _y; }

    ComplexNumber conj() const { return ComplexNumber(_x, -_y); }

    double sqr() const { return _x*_x + _y*_y; }
    double abs() const { return std::sqrt(this->sqr()); }
    double phi() const { return std::atan2(_y, _x); }

    ComplexNumber operator-() const
    {
        return ComplexNumber(-_x, -_y);
    }

    ComplexNumber operator+(const ComplexNumber& z) const
    {
        return ComplexNumber(_x + z.re(), _y + z.im());
    }

    ComplexNumber operator-(const ComplexNumber& z) const
    {
        return ComplexNumber(_x - z.re(), _y - z.im());
    }

    ComplexNumber operator*(const ComplexNumber& z) const
    {
        return ComplexNumber(_x * z.re() - _y * z.im(), _x * z.im() + _y * z.re());
    }

    ComplexNumber operator/(const ComplexNumber& z) const
    {
        double square = z.sqr();
        if (square == 0)
            return ComplexNumber(std::nan("0"), std::nan("0"));

        double real = _x * z.re() + _y * z.im();
        double imag = -_x * z.im() + _y * z.re();

        return ComplexNumber(real / square, imag / square);
    }

    ComplexNumber operator^(const int exponent) const
    {
        double r = std::pow(this->abs(), exponent);
        double f = this->phi() * exponent;

        double real = r * std::cos(f);
        double imag = r * std::sin(f);

        return ComplexNumber(real, imag);
    }

    friend std::ostream& operator<<(std::ostream& os, const ComplexNumber& z)
    {
        os << z.re() << " " << z.im();
        return os;
    }

    bool operator==(const ComplexNumber& z) const
    {
        return (this->re() == z.re()) && (this->im() == z.im());
    }

    bool operator!=(const ComplexNumber& z) const
    {
        return !((*this) == z);
    }

private:
    double _x;
    double _y;
};

using DoublePair = std::pair<double, double>;
using PairVector = std::vector<DoublePair>;

class ComplexVector
{
public:
    friend class FftImplementationUnitTests;
    friend ComplexVector operator/(const ComplexNumber& z, const ComplexVector& vec);

    ComplexVector(const unsigned nsamp)
    {
        _vec.resize(nsamp);
    }

    ComplexVector(const unsigned nsamp, const ComplexNumber& z)
    {
        _vec.resize(nsamp);
        for (unsigned i = 0; i < _vec.size(); i++)
        {
            _vec[i] = z;
        }
    }

    ComplexVector(const unsigned nsamp, const double d)
    {
        _vec.resize(nsamp);
        ComplexNumber z(d, 0.0);
        for (unsigned i = 0; i < _vec.size(); i++)
        {
            _vec[i] = z;
        }
    }

    ComplexVector(PairVector pvec)
    {
        _vec.resize(pvec.size());
        for (unsigned i = 0; i < _vec.size(); i++)
        {
            _vec[i] = ComplexNumber(pvec[i].first, pvec[i].second);
        }
    }

    ComplexVector()
    {
        _vec.resize(0);
    }

    unsigned size() const { return _vec.size(); };

    ComplexNumber operator[](unsigned i) const
    {
        if (i >= this->size())
        {
            std::stringstream ss;
            ss << "Error: ComplexVector access: "
               << i << " >= " << this->size();
            throw std::invalid_argument(ss.str());
        }

        return _vec[i];
    }

private:
    ComplexNumber& operator[](unsigned i)
    {
        if (i >= this->size())
        {
            std::stringstream ss;
            ss << "Error: ComplexVector ref access: "
               << i << " >= " << this->size();
            throw std::invalid_argument(ss.str());
        }

        return _vec[i];
    }

public:
    PairVector getPairVector()
    {
        PairVector pvec(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            pvec[i].first = _vec[i].re();
            pvec[i].second = _vec[i].im();
        }
        return pvec;
    }

    bool setElement(
            const unsigned i,
            const double real,
            const double imag)
    {
        if (i >= this->size())
            return false;

        _vec[i] = ComplexNumber(real, imag);
        return true;
    }

    bool setElement(
            const unsigned i,
            const ComplexNumber z)
    {
        if (i >= this->size())
            return false;

        _vec[i] = z;
        return true;
    }

    bool getElement(
            const unsigned i,
            ComplexNumber& z) const
    {
        if (i >= this->size())
            return false;

        z = _vec[i];
        return true;
    }

    std::vector<double> re() const
    {
        std::vector<double> vec;
        for (const ComplexNumber& z : _vec)
        {
            vec.emplace_back(z.re());
        }
        return vec;
    }

    std::vector<double> im() const
    {
        std::vector<double> vec;
        for (const ComplexNumber& z : _vec)
        {
            vec.emplace_back(z.im());
        }
        return vec;
    }

    ComplexVector conj() const
    {
        ComplexVector vec(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            vec[i] = _vec[i].conj();
        }
        return vec;
    }

    std::vector<double> sqr() const
    {
        std::vector<double> out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i].sqr();
        }
        return out;
    }
    
    std::vector<double> abs() const
    {
        std::vector<double> out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i].abs();
        }
        return out;
    }
    
    std::vector<double> phi() const
    {
        std::vector<double> out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i].phi();
        }
        return out;
    }
    
    ComplexVector operator-() const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = -_vec[i];
        }
        return out;
    }

    ComplexVector operator+(const ComplexVector& vec) const
    {
        if (vec.size() != _vec.size())
            return ComplexVector(0);

        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] + vec[i];
        }
        return out;
    }

    ComplexVector operator-(const ComplexVector& vec) const
    {
        if (vec.size() != _vec.size())
            return ComplexVector(0);

        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] - vec[i];
        }
        return out;
    }

    ComplexVector operator*(const ComplexVector& vec) const
    {
        if (vec.size() != _vec.size())
            return ComplexVector(0);

        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] * vec[i];
        }
        return out;
    }

    ComplexVector operator/(const ComplexVector& vec) const
    {
        if (vec.size() != _vec.size())
            return ComplexVector(0);

        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] / vec[i];
        }
        return out;
    }

    ComplexVector operator+(const ComplexNumber& z) const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] + z;
        }
        return out;
    }

    ComplexVector operator-(const ComplexNumber& z) const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] - z;
        }
        return out;
    }

    ComplexVector operator*(const ComplexNumber& z) const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] * z;
        }
        return out;
    }

    ComplexVector operator/(const ComplexNumber& z) const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i] / z;
        }
        return out;
    }

    ComplexVector operator^(const int exponent) const
    {
        ComplexVector out(this->size());
        for (unsigned i = 0; i < this->size(); i++)
        {
            out[i] = _vec[i]^exponent;
        }
        return out;
    }

    bool operator==(const ComplexVector& vec) const
    {
        for (unsigned i = 0; i < this->size(); i++)
        {
            if (_vec[i] != vec[i])
                return false;
        }
        return true;
    }

    bool operator!=(const ComplexVector& vec) const
    {
        return !((*this) == vec);
    }

    bool operator==(const ComplexNumber& z) const
    {
        for (unsigned i = 0; i < this->size(); i++)
        {
            if (_vec[i] != z)
                return false;
        }
        return true;
    }

    bool operator!=(const ComplexNumber& z) const
    {
        return !((*this) == z);
    }

    friend std::ostream& operator<<(std::ostream& os, const ComplexVector& vec)
    {
        for (unsigned i = 0; i < vec.size(); i++)
            os << vec[i] << std::endl;
        return os;
    }

private:
    std::vector<ComplexNumber> _vec;
};

bool operator==(const ComplexNumber& z, const double d)
{
    return z == ComplexNumber(d, 0.0);
}

bool operator==(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) == z;
}

bool operator!=(const ComplexNumber& z, const double d)
{
    return z != ComplexNumber(d, 0.0);
}

bool operator!=(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) != z;
}

bool operator==(const ComplexNumber& z, const ComplexVector& vec)
{
    return vec == z;
}

bool operator!=(const ComplexNumber& z, const ComplexVector& vec)
{
    return vec != z;
}

bool operator==(const ComplexVector& vec, const double d)
{
    return vec == ComplexNumber(d, 0.0);
}

bool operator==(const double d, const ComplexVector& vec)
{
    return vec == ComplexNumber(d, 0.0);
}

bool operator!=(const ComplexVector& vec, const double d)
{
    return vec != ComplexNumber(d, 0.0);
}

bool operator!=(const double d, const ComplexVector& vec)
{
    return vec != ComplexNumber(d, 0.0);
}

ComplexVector operator+(const ComplexNumber& z, const ComplexVector& vec)
{
    return vec + z;
}

ComplexVector operator-(const ComplexNumber& z, const ComplexVector& vec)
{
    return -vec + z;
}

ComplexVector operator*(const ComplexNumber& z, const ComplexVector& vec)
{
    return vec * z;
}

ComplexVector operator/(const ComplexNumber& z, const ComplexVector& vec)
{
        ComplexVector out(vec.size());
        for (unsigned i = 0; i < vec.size(); i++)
        {
            out[i] = z / vec[i];
        }
        return out;    
}

ComplexNumber operator+(const ComplexNumber& z, const double d)
{
    return z + ComplexNumber(d, 0.0);
}

ComplexNumber operator+(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) + z;
}

ComplexNumber operator-(const ComplexNumber& z, const double d)
{
    return z - ComplexNumber(d, 0.0);
}

ComplexNumber operator-(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) - z;
}

ComplexNumber operator*(const ComplexNumber& z, const double d)
{
    return z * ComplexNumber(d, 0.0);
}

ComplexNumber operator*(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) * z;
}

ComplexNumber operator/(const ComplexNumber& z, const double d)
{
    return z / ComplexNumber(d, 0.0);
}

ComplexNumber operator/(const double d, const ComplexNumber& z)
{
    return ComplexNumber(d, 0.0) / z;
}

ComplexVector operator+(const ComplexVector& vec, const double d)
{
    return vec + ComplexNumber(d, 0.0);
}

ComplexVector operator+(const double d, const ComplexVector& vec)
{
    return ComplexNumber(d, 0.0) + vec;
}

ComplexVector operator-(const ComplexVector& vec, const double d)
{
    return vec - ComplexNumber(d, 0.0);
}

ComplexVector operator-(const double d, const ComplexVector& vec)
{
    return ComplexNumber(d, 0.0) - vec;
}

ComplexVector operator*(const ComplexVector& vec, const double d)
{
    return vec * ComplexNumber(d, 0.0);
}

ComplexVector operator*(const double d, const ComplexVector& vec)
{
    return ComplexNumber(d, 0.0) * vec;
}

ComplexVector operator/(const ComplexVector& vec, const double d)
{
    return vec / ComplexNumber(d, 0.0);
}

ComplexVector operator/(const double d, const ComplexVector& vec)
{
    return ComplexNumber(d, 0.0) / vec;
}

#endif
