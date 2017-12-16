#include "sweepline.hpp"

#include <utility>
#include <limits>
#include <iterator>
#include <algorithm>
#include <random>
#include <tuple>
#include <set>
#include <vector>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#ifndef NDEBUG
#include <unordered_set>
#endif

#include <cassert>
#include <cmath>

template< typename iterator >
struct proxy_iterator
{

    using iterator_type     = typename std::iterator_traits< iterator >::value_type;
    using iterator_traits   = std::iterator_traits< iterator_type >;

    using iterator_category = std::forward_iterator_tag;
    using value_type        = typename iterator_traits::value_type;
    using difference_type   = typename iterator_traits::difference_type;
    using pointer           = typename iterator_traits::pointer;
    using reference         = typename iterator_traits::reference;

    iterator it;

    operator iterator_type () const { return *it; }

    proxy_iterator & operator ++ () { ++it; return *this; }
    proxy_iterator operator ++ (int) { return {it++}; }

    reference operator * () const { return **it; }
    pointer operator -> () const { return &operator * (); }

    bool operator == (const proxy_iterator pi) const { return (it == pi.it); }
    bool operator != (const proxy_iterator pi) const { return (it != pi.it); }

};

template< typename point, typename value_type = decltype(std::declval< point >().x) >
struct voronoi
{

    using size_type = std::size_t;

    bool draw_indices = false;
    bool draw_circles = false;

    value_type eps = value_type(10) * std::numeric_limits< value_type >::epsilon();
    value_type delta = value_type(0.001);

    std::ostream & log_;

    explicit
    voronoi(std::ostream & _log)
        : log_(_log)
    {
        assert(!(delta < eps));
        log_ << "eps = " << eps << '\n';
        log_ << "delta = " << delta << '\n';
    }

private :

    const value_type zero = value_type(0);
    const value_type one = value_type(1);

    std::mt19937 rng;
    std::normal_distribution< value_type > normal_;
    std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};

public :

    using seed_type = typename std::mt19937::result_type;

    void seed(const seed_type seed)
    {
        log_ << "seed = " << seed << '\n';
        rng.seed(seed);
    }

    struct less
    {

        const value_type & eps;

        bool operator () (const value_type & l,
                          const value_type & r) const
        {
            return l + eps < r;
        }

        bool operator () (const value_type & lx, const value_type & ly,
                          const value_type & rx, const value_type & ry) const
        {
            if (operator () (lx, rx)) {
                return true;
            } else if (operator () (rx, lx)) {
                return false;
            } else {
                return operator () (ly, ry);
            }
        }

        bool operator () (const point & l, const point & r) const
        {
            return operator () (l.x, l.y, r.x, r.y);
        }

    };

    void ball(std::ostream & _out, const value_type radius, const size_type N)
    {
        std::set< point, less > unique_points_{less{delta}};
        _out << N << '\n';
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
            size_type m = 0;
            do {
                point p{normal_(rng), normal_(rng)};
                value_type norm = p.x * p.x + p.y * p.y;
                if (eps * eps < norm) {
                    using std::sqrt;
                    norm = radius * sqrt(zero_to_one_(rng) / std::move(norm));
                    p.x *= norm;
                    p.y *= norm;
                } else {
                    p.x = p.y = zero;
                }
                if (unique_points_.insert(std::move(p)).second) {
                    break;
                }
            } while (++m < M);
            if (m == M) {
                log_ << "the number (" << M << ") of attempts is exceeded\n";
                log_ << "only " << n << "points generated\n";
                break;
            }
        }
        for (const point & point_ : unique_points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    void square(std::ostream & _out, const value_type bbox, const size_type N)
    {
        std::set< point, less > unique_points_{less{delta}};
        _out << N << '\n';
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed square
            size_type m = 0;
            do {
                point p{zero_to_one_(rng), zero_to_one_(rng)};
                p.x += p.x;
                p.y += p.y;
                p.x -= one;
                p.y -= one;
                p.x *= bbox;
                p.y *= bbox;
                if (unique_points_.insert(std::move(p)).second) {
                    break;
                }
            } while (++m < M);
            if (m == M) {
                log_ << "the number (" << M << ") of attempts is exceeded\n";
                log_ << "only " << n << "points generated\n";
                break;
            }
        }
        for (const point & point_ : unique_points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    void rectangular_grid(std::ostream & _out, const size_type bbox) const
    {
        const size_type N = 1 + 4 * bbox * (bbox + 1);
        _out << N << '\n';
        _out << "0 0\n";
        for (size_type x = 1; x <= bbox; ++x) {
            _out << "0 " << x << '\n';
            _out << x << " 0\n";
            _out << "0 -" << x << '\n';
            _out << '-' << x << " 0\n";
            for (size_type y = 1; y <= bbox; ++y) {
                _out << x << ' ' << y << '\n';
                _out << '-' << x << ' ' << y << '\n';
                _out << x << " -" << y << '\n';
                _out << '-' << x << " -" << y << '\n';
            }
        }
    }

    void diagonal_grid(std::ostream & _out, const size_type bbox) const
    {
        const size_type N = (1 + 2 * bbox * (bbox + 1));
        _out << N << '\n';
        _out << "0 0\n";
        size_type i = 1;
        for (size_type x = 1; x <= bbox; ++x) {
            for (size_type y = x % 2; y <= bbox; y += 2) {
                _out << x << ' ' << y << '\n';
                _out << y << " -" << x << '\n';
                _out << '-' << x << " -" << y << '\n';
                _out << '-' << y << ' ' << x << '\n';
                i += 4;
            }
        }
        assert(N == i);
    }

    void hexagonal_grid(std::ostream & _out, const size_type size) const
    {
        const size_type N = (size + size) * (size + 1);
        _out << N << '\n';
        using std::sqrt;
        const value_type step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            const value_type xx = value_type(x) * step;
            if ((x % 2) == 0) {
                const value_type yy = value_type(3 * (x - 1));
                _out << "0 " << yy << '\n';
                _out << "0 -" << yy << '\n';
            } else {
                _out << xx << " 0\n";
                _out << '-' << xx << " 0\n";
            }
            i += 2;
            for (size_type y = 1 + (x % 2); y <= size; y += 2) {
                const size_type yy = 3 * y;
                _out << xx << ' ' << yy << '\n';
                _out << xx << " -" << yy << '\n';
                _out << '-' << xx << ' ' << yy << '\n';
                _out << '-' << xx << " -" << yy << '\n';
                i += 4;
            }
        }
        if ((size % 2) != 0) {
            const value_type yy = value_type(3 * size);
            _out << "0 " << yy << '\n';
            _out << "0 -" << yy << '\n';
            i += 2;
        }
        assert(i == N);
    }

    void triangular_grid(std::ostream & _out, const size_type size) const
    {
        const size_type N = (1 + size + size) * (size + size);
        _out << N << '\n';
        using std::sqrt;
        const value_type step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            const size_type xx = 3 * x;
            {
                const value_type yy = value_type(step) * (xx - 1 - (x % 2));
                _out << "0 " << yy << '\n';
                _out << "0 -" << yy << '\n';
            }
            i += 2;
            for (size_type y = 1; y <= size; ++y) {
                const value_type yy = step * value_type(3 * y - 1 - ((y + x) % 2));
                _out << xx << ' ' << yy << '\n';
                _out << xx << " -" << yy << '\n';
                _out << '-' << xx << ' ' << yy << '\n';
                _out << '-' << xx << " -" << yy << '\n';
                i += 4;
            }
        }
        assert(i == N);
    }

    struct ipoint { size_type x, y; };

    template< std::size_t nsqr >
    void quadrant(std::ostream & _out, const size_type max, const ipoint (& q)[nsqr]) const
    {
        if (0 == max) {
            _out << (nsqr * 8) << '\n';
        } else {
            _out << (nsqr * 8 + 4) << '\n';
            _out << "0 " << max << '\n';
            _out << "0 -" << max << '\n';
            _out << max << " 0\n";
            _out << '-' << max << " 0\n";
        }
        const auto qprint = [&] (const bool swap, const bool sx, const bool sy)
        {
            for (const ipoint & p : q) {
                assert(p.x != 0);
                assert(p.y != 0);
                if (sx) {
                    _out << '-';
                }
                _out << (swap ? p.x : p.y) << ' ';
                if (sy) {
                    _out << '-';
                }
                _out << (swap ? p.y : p.x) << '\n';
            }
        };
        qprint(true,  true,  true);
        qprint(true,  false, true);
        qprint(true,  true,  false);
        qprint(true,  false, false);
        qprint(false, true,  true);
        qprint(false, false, true);
        qprint(false, true,  false);
        qprint(false, false, false);
    }

    using points = std::vector< point >;

private :

    points points_;

    void input(std::istream & _in)
    {
        size_type M = 0;
        if (!(_in >> M)) {
            assert(false);
        }
        points_.reserve(M);
        for (size_type m = 0; m < M; ++m) {
            points_.emplace_back();
            point & point_ = points_.back();
            if (!(_in >> point_.x)) {
                assert(false);
            }
            if (!(_in >> point_.y)) {
                assert(false);
            }
        }
    }

public :

    friend
    std::istream & operator >> (std::istream & _in, voronoi & _voronoi)
    {
        _voronoi.input(_in);
        return _in;
    }

    void swap_xy()
    {
        for (point & point_ : points_) {
            using std::swap;
            swap(point_.x, point_.y);
        }
    }

    void shift_xy(const value_type & dx, const value_type & dy)
    {
        for (point & point_ : points_) {
            point_.x += dx;
            point_.y += dy;
        }
    }

    void rotate(const value_type & angle)
    {
        using std::cos;
        using std::sin;
        const value_type cosine = cos(angle);
        const value_type sine = sin(angle);
        for (point & point_ : points_) {
            value_type x = cosine * point_.x - sine * point_.y;
            point_.y = sine * point_.x + cosine * point_.y;
            point_.x = std::move(x);
        }
    }

    using site = typename points::const_iterator;

    using sweepline_type = sweepline< site, point, value_type >;

    sweepline_type sweepline_{eps};

    void operator () ()
    {
        assert((std::set< point, less >{std::cbegin(points_), std::cend(points_), less{delta}}.size() == points_.size()));
        log_ << "N = " << points_.size() << '\n';
#if 0
        std::sort(std::begin(points_), std::end(points_));
        sweepline_(std::cbegin(points_), std::cend(points_));
#else
        using sites = std::vector< site >;
        sites sites_;
        {
            sites_.reserve(points_.size() + 1);
            const auto send = std::cend(points_);
            for (auto s = std::cbegin(points_); s != send; ++s) {
                sites_.push_back(s);
            }
            const auto sless = [] (const site l, const site r) -> bool { return *l < *r; };
            std::sort(std::begin(sites_), std::end(sites_), sless);
            sites_.push_back(send);
        }
        using psite = proxy_iterator< typename sites::const_iterator >;
        sweepline_(psite{std::cbegin(sites_)}, psite{std::prev(std::cend(sites_))});
#ifndef NDEBUG
        using pvertex = typename sweepline_type::pvertex;
        const auto vhash = [] (const pvertex v)
        {
            return std::hash< typename std::iterator_traits< pvertex >::pointer >{}(&*v);
        };
        using vpoints_type = std::unordered_multiset< pvertex, decltype((vhash)) >;
        const size_type capacity_ = 1 + (sweepline_.vertices_.size() * 3) / 2;
        vpoints_type heads{capacity_, vhash};
        vpoints_type tails{capacity_, vhash};
        for (const auto & edge_ : sweepline_.edges_) {
            if (edge_.b != sweepline_.nv) {
                heads.insert(edge_.b);
            }
            if (edge_.e != sweepline_.nv) {
                tails.insert(edge_.e);
            }
        }
        for (auto v = std::begin(sweepline_.vertices_); v != std::end(sweepline_.vertices_); ++v) {
            const size_type b = heads.count(v);
            assert(0 < b);
            assert(2 < tails.count(v) + b);
        }
#endif
#endif
    }

    struct truncate_edge
    {

        const point & l;
        const point & r;
        const point & p;

        const point & vmin;
        const point & vmax;

        const value_type & eps_;

        const value_type dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
        const value_type dy = l.x - r.x;

        point px(const value_type & y) const { return {(p.x + (y - p.y) * dx / dy), y}; }

        point py(const value_type & x) const
        {
            const value_type y = p.y + (x - p.x) * dy / dx;
            if (+eps_ < dy) {
                if (vmax.y < y) {
                    return px(vmax.y);
                }
            } else if (dy < -eps_) {
                if (y < vmin.y) {
                    return px(vmin.y);
                }
            }
            return {x, y};
        }

        operator point () const
        {
            if (+eps_ < dx) {
                return py(vmax.x);
            } else if (dx < -eps_) {
                return py(vmin.x);
            } else {
                if (+eps_ < dy) {
                    return {p.x, vmax.y};
                } else if (dy < -eps_) {
                    return {p.x, vmin.y};
                } else {
                    assert(false);
                    return {p.x, p.y};
                }
            }
        }

    };

    value_type zoom = value_type(0.2);

    template< typename V, typename E >
    void output(std::ostream & _gnuplot,
                const V & _vertices,
                const E & _edges) const
    {
        if (points_.empty()) {
            _gnuplot << "print 'no point to process'\n;";
            return;
        }
        // pmin, pmax denotes bounding box
        point vmin = points_.front();
        point vmax = vmin;
        const auto pminmax = [&] (const point & p)
        {
            if (p.x < vmin.x) {
                vmin.x = p.x;
            } else if (vmax.x < p.x) {
                vmax.x = p.x;
            }
            if (p.y < vmin.y) {
                vmin.y = p.y;
            } else if (vmax.y < p.y) {
                vmax.y = p.y;
            }
        };
        std::for_each(std::next(std::cbegin(points_)), std::cend(points_), pminmax);
        assert(value_type(-0.5) < zoom);
        if (vmin.x + delta < vmax.x) {
            value_type dx = vmax.x - vmin.x;
            dx *= zoom;
            vmin.x -= dx;
            vmax.x += dx;
        } else {
            vmin.x -= delta;
            vmax.x += delta;
        }
        if (vmin.y + delta < vmax.y) {
            value_type dy = vmax.y - vmin.y;
            dy *= zoom;
            vmin.y -= dy;
            vmax.y += dy;
        } else {
            vmin.y -= delta;
            vmax.y += delta;
        }
        point pmin = vmin;
        point pmax = vmax;
        std::for_each(std::cbegin(_vertices), std::cend(_vertices), [&] (const auto & v) { pminmax(v.c); });
        {
            _gnuplot << "set size square;\n"
                        "set key left;\n"
                        "unset colorbox;\n"
                        "set cbrange [-1:1];\n";
            _gnuplot << "set xrange [" << pmin.x << ':' << pmax.x << "];\n";
            _gnuplot << "set yrange [" << pmin.y << ':' << pmax.y << "];\n";
            _gnuplot << "set size ratio -1;\n";
        }
        {
            _gnuplot << "$sites << EOI\n";
            size_type i = 0;
            for (const point & point_ : points_) {
                _gnuplot << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
            }
            _gnuplot << "EOI\n";
        }
        if (draw_circles && !_vertices.empty()) {
            _gnuplot << "$circles << EOI\n";
            size_type i = 0;
            for (const auto & vertex_ : _vertices) {
                _gnuplot << vertex_.c.x << ' ' << vertex_.c.y << ' ' << vertex_.R << ' ' << i++ << '\n';
            }
            _gnuplot << "EOI\n";
        }
        if (!_edges.empty()) {
            const auto pout = [&] (const point & p)
            {
                _gnuplot << p.x << ' ' << p.y << '\n';
            };
            _gnuplot << "$edges << EOI\n";
            const auto nv = std::end(_vertices);
            for (const auto & edge_ : _edges) {
                const bool beg = (edge_.b != nv);
                const bool end = (edge_.e != nv);
                const point & l = *edge_.l;
                const point & r = *edge_.r;
                if (beg != end) {
                    const point & p = (beg ? edge_.b : edge_.e)->c;
                    if (!(p.x < vmin.x) && !(vmax.x < p.x) && !(p.y < vmin.y) && !(vmax.y < p.y)) {
                        pout(p);
                        pout(truncate_edge{(beg ? l : r), (end ? l : r), p, vmin, vmax, eps});
                    }
                } else if (beg) {
                    assert(!less{eps}(edge_.e->c, edge_.b->c));
                    pout(edge_.b->c);
                    pout(edge_.e->c);
                } else {
                    const point p{(l.x + r.x) / value_type(2), (l.y + r.y) / value_type(2)};
                    pout(truncate_edge{l, r, p, vmin, vmax, eps});
                    pout(truncate_edge{r, l, p, vmin, vmax, eps});
                }
                _gnuplot << "\n"; // separate lines
            }
            _gnuplot << "EOI\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '$sites' with points title 'sites # " << points_.size() << "'";
        if (draw_indices) {
            _gnuplot << ", '' with labels offset character 0, character 1 notitle";
        }
        if (draw_circles && !_vertices.empty()) {
            _gnuplot << ", '$circles' with circles title 'vertices # " << _vertices.size() << "' linecolor palette";
        }
        if (!_edges.empty()) {
            _gnuplot << ", '$edges' with lines title 'edges # " << _edges.size() <<  "'";
        }
        _gnuplot << ";\n";
    }

private :

    void output(std::ostream & _out) const
    {
        return output(_out, sweepline_.vertices_, sweepline_.edges_);
    }

public :

    friend
    std::ostream & operator << (std::ostream & _gnuplot, const voronoi & _voronoi)
    {
        _voronoi.output(_gnuplot);
        return _gnuplot;
    }

};

#include <iterator>
#include <algorithm>
#include <limits>
#include <vector>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <chrono>
#include <typeinfo>
#include <memory>
#include <mutex>

#include <cstdlib>

#include <x86intrin.h>
#include <cxxabi.h>

#define RED(str) __extension__ "\e[1;31m" str "\e[0m"

namespace
{

#pragma GCC diagnostic push
#ifdef __clang__
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif
std::mutex m;
std::unique_ptr< char, decltype((std::free)) > demangled_name{nullptr, std::free};
std::size_t length = 0;
#pragma GCC diagnostic pop

}

inline
std::string
get_demangled_name(const char * const symbol) noexcept
{
    if (!symbol) {
        return "<null>";
    }
    std::lock_guard< std::mutex > lock(m);
    int status = -4;
    demangled_name.reset(abi::__cxa_demangle(symbol, demangled_name.release(), &length, &status));
    return ((status == 0) ? demangled_name.get() : symbol);
}

using value_type = double;

struct alignas(__m128d) point
{

    value_type x, y;

    bool operator < (const point & p) const
    {
        return std::tie(x, y) < std::tie(p.x, p.y);
    }

};

inline
std::ostream & operator << (std::ostream & out, const point & p)
{
    return out << p.x << ' ' << p.y;
}

int main()
{
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    using std::chrono::steady_clock;

    using voronoi_type = voronoi< point >;
    std::ostream & log_ = std::clog;
    voronoi_type voronoi_{log_};
    { // input
#if 0
        std::istream & in_ = std::cin;
#else
        std::stringstream in_;
# ifdef _LIBCPP_VERSION
        in_ >> std::hexfloat;
# else
        in_.precision(std::numeric_limits< value_type >::digits10 + 1);
# endif
# if 0
        voronoi_.draw_circles = true;
        voronoi_.draw_indices = true;
#  if 0
        in_ << "0\n";
#  elif 0
        in_ << "1\n"
               "0 0\n";
#  elif 0
        in_ << "2\n"
               "0 0\n"
               "1 0\n";
#  elif 0
        in_ << "2\n"
               "0 0\n"
               "0 1\n";
#  elif 0
        in_ << "2\n"
               "0 0\n"
               "1 1\n";
#  elif 0
        in_ << "2\n"
               "0 0\n"
               "1 -1\n";
#  elif 0
        in_ << "3\n"
               "0 0\n"
               "0 1\n"
               "0 2\n";
#  elif 0
        in_ << "5\n"
               "0 0\n"
               "0 1\n"
               "0 2\n"
               "0 3\n"
               "0 4\n";
#  elif 0
        in_ << "5\n"
               "0 0\n"
               "1 0\n"
               "2 0\n"
               "3 0\n"
               "4 0\n";
#  elif 0
        in_ << "5\n"
               "0 0\n"
               "1 1\n"
               "2 2\n"
               "3 3\n"
               "4 4\n";
#  elif 0
        in_ << "5\n"
               "0 -0\n"
               "1 -1\n"
               "2 -2\n"
               "3 -3\n"
               "4 -4\n";
#  elif 0
        in_ << "3\n"
               "0 0\n"
               "1 -1\n"
               "1 1\n";
#  elif 0
        in_ << "3\n"
               "-1 0\n"
               "0 1\n"
               "1 0\n";
#  elif 0
        in_ << "3\n"
               "-1 0\n"
               "0 -1\n"
               "1 0\n";
#  elif 0
        in_ << "3\n"
               "0 -1\n"
               "0 1\n"
               "1 0\n";
#  elif 0
        in_ << "4\n"
               "-3 -4\n"
               "-3 4\n"
               "-4 -3\n"
               "-4 3\n";
#  else
        // very good test!
        in_ << "4\n"
               "1 2\n"
               "2 3\n"
               "3 0\n"
               "3 2\n"
               ;
#  endif
# elif 0
        // Concentric:
        voronoi_.draw_circles = true;
#  if 0
        in_ << "4\n"
               "-1 0\n"
               "0 -1\n"
               "0 1\n"
               "1 0\n";
#  elif 0
        voronoi_.quadrant(in_, /*5*/0, {{3, 4}});
#  elif 0
        voronoi_.quadrant(in_, 25, {{7, 24}, {15, 20}});
#  elif 0
        voronoi_.quadrant(in_, 65, {{16, 63}, {25, 60}, {33, 56}, {39, 52}});
#  elif 0
        voronoi_.quadrant(in_, 325, {{36, 323}, {80, 315}, {91, 312}, {125, 300}, {165, 280}, {195, 260}, {204, 253}});
#  elif 0
        voronoi_.quadrant(in_, 1105, {{47,  1104}, {105, 1100}, {169, 1092}, {264, 1073}, {272, 1071}, {425, 1020}, {468, 1001},
                                      {520, 975 }, {561, 952 }, {576, 943 }, {663, 884 }, {700, 855 }, {744, 817 }});
#  else
        voronoi_.quadrant(in_, 5525, {{235,  5520}, {525,  5500}, {612,  5491}, {845,  5460},
                                      {1036, 5427}, {1131, 5408}, {1320, 5365}, {1360, 5355}, {1547, 5304},
                                      {2044, 5133}, {2125, 5100}, {2163, 5084}, {2340, 5005}, {2600, 4875},
                                      {2805, 4760}, {2880, 4715}, {3124, 4557}, {3315, 4420}, {3468, 4301},
                                      {3500, 4275}, {3720, 4085}, {3861, 3952}});
#  endif
# else
        // Rectangular, rectangular pi / 4 rotated, hexagonal, triangular grids and points uniformely distributed into a circle or square:
        {
            using seed_type = typename voronoi_type::seed_type;
#  if 0
            const seed_type seed = 2911579113;
#  else
            std::random_device D;
            const auto seed = static_cast< seed_type >(D());
#  endif
            voronoi_.seed(seed);
            //gnuplot_ << "set title 'seed = 0x" << std::hex << std::nouppercase << seed << ", N = " <<  std::dec << N << "'\n";
        }
#  if 0
        voronoi_.rectangular_grid(in_, 10); voronoi_.draw_circles = true;
#  elif 0
        voronoi_.diagonal_grid(in_, 20); voronoi_.draw_circles = true;
#  elif 0
        voronoi_.hexagonal_grid(in_, 20); //voronoi_.eps = value_type(0.0001); //voronoi_.draw_circles = true;
#  elif 0
        voronoi_.triangular_grid(in_, 24); voronoi_.eps = value_type(0.000000001); //voronoi_.draw_circles = true;
#  elif 0
        voronoi_.square(in_, value_type(10000), 100000);
#  else
        voronoi_.ball(in_, value_type(10000), 100000); // voronoi_.draw_circles = true; // voronoi_.draw_indices = true;
#  endif
# endif
        //log_ << in_.str() << '\n';
#endif
        {
            const auto start = steady_clock::now();
            in_ >> voronoi_;
            log_ << "input time = "
                 << duration_cast< microseconds >(steady_clock::now() - start).count()
                 << "us\n";
        }
        //voronoi_.rotate(M_PI_4);
        //voronoi_.swap_xy();
        //voronoi_.shift_xy(value_type(10000), value_type(10000));
    }
    { // run
        log_ << "start\n";
        const auto start = steady_clock::now();
        try {
            for (std::size_t i = 0; i < 1; ++i) {
                voronoi_.sweepline_.clear();
                voronoi_();
            }
        } catch (...) {
            log_ <<  RED("Exception catched!") "\n";
            if (std::type_info * et = abi::__cxa_current_exception_type()) {
                log_ << "exception type: " << get_demangled_name(et->name()) << '\n';
            } else {
                log_ << "unknown exception\n";
            }
        }
        log_ << "sweepline time = "
             << duration_cast< microseconds >(steady_clock::now() - start).count()
             << "us\n";
    }
    { // output
        //voronoi_.draw_circles = false; // (sweepline_.vertices_.size() < 300);
        //voronoi_.draw_indices = false;
        using sweepline_type = typename voronoi_type::sweepline_type;
        const sweepline_type & sweepline_ = voronoi_.sweepline_;
        log_ << "vertices # " << sweepline_.vertices_.size() << '\n';
        log_ << "edges # " << sweepline_.edges_.size() << '\n';
        std::string command_line_ = "sweepline.plt";
        std::ofstream f{command_line_};
        std::ostream & gnuplot_ = f;//std::cout;
#if 0
        { // clone (O(|vertices| * |edges|))
            using vertices = std::vector< typename sweepline_type::vertex >;
            using pvertex = typename vertices::const_iterator;
            using site = typename voronoi_type::site;
            const vertices vertices_{std::cbegin(sweepline_.vertices_), std::cend(sweepline_.vertices_)};
            const pvertex nv = std::cend(vertices_);
            const auto vclone = [&] (typename sweepline_type::const pvertex v) -> pvertex
            {
                return std::prev(nv, std::distance(v, sweepline_.nv));
            };
            struct edge { site l, r; pvertex b, e; };
            const auto eclone = [&] (typename sweepline_type::const edge & _edge) -> edge
            {
                return {_edge.l, _edge.r, vclone(_edge.b), vclone(_edge.e)};
            };
            std::vector< edge > edges_;
            edges_.reserve(sweepline_.edges_.size());
            std::transform(std::cbegin(sweepline_.edges_), std::cend(sweepline_.edges_),
                           std::back_inserter(edges_),
                           eclone);
            voronoi_.output(gnuplot_, vertices_, edges_);
            gnuplot_ << std::endl;
        }
#else
        gnuplot_ << voronoi_ << std::endl;
#endif
        command_line_.insert(0, "gnuplot -p ");
#ifndef _WIN32
#if 0
        command_line_.insert(0, "GNUTERM=qt ");
#elif 0
        command_line_.insert(0, "GNUTERM=wxt ");
#endif
#endif
        log_ << std::flush;
        gnuplot_ << std::flush;
        if (std::system(command_line_.c_str()) != 0) {
            log_ << "error occured during execution of \"" << command_line_ << "\" command line\n";
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
