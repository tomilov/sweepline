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

    proxy_iterator(iterator const i) : it{i} { ; }
    operator iterator_type () const { return *it; }

    proxy_iterator & operator ++ () { ++it; return *this; }
    proxy_iterator operator ++ (int) { return {it++}; }

    reference operator * () const { return **it; }
    pointer operator -> () const { return &operator * (); }

    bool operator == (proxy_iterator const & i) const { return (it == i.it); }
    bool operator != (proxy_iterator const & i) const { return !operator == (i); }

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

    value_type const zero = value_type(0);
    value_type const one = value_type(1);

    std::mt19937 rng;
    std::normal_distribution< value_type > normal_;
    std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};

public :

    using seed_type = typename std::mt19937::result_type;

    void
    seed(seed_type const seed)
    {
        log_ << "seed = " << seed << '\n';
        rng.seed(seed);
    }

    struct less
    {

        value_type const & eps;

        bool operator () (value_type const & l,
                          value_type const & r) const
        {
            return l + eps < r;
        }

        bool operator () (value_type const & lx, value_type const & ly,
                          value_type const & rx, value_type const & ry) const
        {
            if (operator () (lx, rx)) {
                return true;
            } else if (operator () (rx, lx)) {
                return false;
            } else {
                return operator () (ly, ry);
            }
        }

        bool operator () (point const & l, point const & r) const
        {
            return operator () (l.x, l.y, r.x, r.y);
        }

    };

    void
    ball(std::ostream & _out, value_type const radius, size_type const N)
    {
        std::set< point, less > unique_points_{less{delta}};
        _out << N << '\n';
        value_type const twosqreps = eps;
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
            size_type m = 0;
            do {
                point p{normal_(rng), normal_(rng)};
                value_type norm = p.x * p.x + p.y * p.y;
                if (twosqreps < norm) {
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
        for (point const & point_ : unique_points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    void
    square(std::ostream & _out, value_type const bbox, size_type const N)
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
        for (point const & point_ : unique_points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    void
    rectangle_grid(std::ostream & _out, size_type const bbox) const
    {
        size_type const N = 1 + 4 * bbox * (bbox + 1);
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

    void
    diagonal_grid(std::ostream & _out, size_type const bbox) const
    {
        size_type const N = (1 + 2 * bbox * (bbox + 1));
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

    void
    hexagonal_grid(std::ostream & _out, size_type const size) const
    {
        size_type const N = (size + size) * (size + 1);
        _out << N << '\n';
        using std::sqrt;
        value_type const step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            value_type const xx = value_type(x) * step;
            if ((x % 2) == 0) {
                value_type const yy = value_type(3 * (x - 1));
                _out << "0 " << yy << '\n';
                _out << "0 -" << yy << '\n';
            } else {
                _out << xx << " 0\n";
                _out << '-' << xx << " 0\n";
            }
            i += 2;
            for (size_type y = 1 + (x % 2); y <= size; y += 2) {
                size_type const yy = 3 * y;
                _out << xx << ' ' << yy << '\n';
                _out << xx << " -" << yy << '\n';
                _out << '-' << xx << ' ' << yy << '\n';
                _out << '-' << xx << " -" << yy << '\n';
                i += 4;
            }
        }
        if ((size % 2) != 0) {
            value_type const yy = value_type(3 * size);
            _out << "0 " << yy << '\n';
            _out << "0 -" << yy << '\n';
            i += 2;
        }
        assert(i == N);
    }

    void
    triangular_grid(std::ostream & _out, size_type const size) const
    {
        size_type const N = (1 + size + size) * (size + size);
        _out << N << '\n';
        using std::sqrt;
        value_type const step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            size_type const xx = 3 * x;
            {
                value_type const yy = value_type(step) * (xx - 1 - (x % 2));
                _out << "0 " << yy << '\n';
                _out << "0 -" << yy << '\n';
            }
            i += 2;
            for (size_type y = 1; y <= size; ++y) {
                value_type const yy = step * value_type(3 * y - 1 - ((y + x) % 2));
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
    void
    quadrant(std::ostream & _out, size_type const max, ipoint const (& q)[nsqr]) const
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
        auto const qprint = [&] (bool const swap, bool const sx, bool const sy)
        {
            for (ipoint const & p : q) {
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
    std::istream &
    operator >> (std::istream & _in, voronoi & _voronoi)
    {
        _voronoi.input(_in);
        return _in;
    }

    void
    swap_xy()
    {
        for (point & point_ : points_) {
            using std::swap;
            swap(point_.x, point_.y);
        }
    }

    void
    shift_xy(value_type const & dx, value_type const & dy)
    {
        for (point & point_ : points_) {
            point_.x += dx;
            point_.y += dy;
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
            auto const send = std::cend(points_);
            for (auto s = std::cbegin(points_); s != send; ++s) {
                sites_.push_back(s);
            }
            auto const sless = [] (site const l, site const r) -> bool { return *l < *r; };
            std::sort(std::begin(sites_), std::end(sites_), sless);
            sites_.push_back(send);
        }
        using psite = proxy_iterator< typename sites::const_iterator >;
        sweepline_(psite{std::cbegin(sites_)}, psite{std::prev(std::cend(sites_))});
#endif
    }

    value_type zoom = value_type(0.2);

    template< typename V, typename E >
    void output(std::ostream & _gnuplot,
                V const & _vertices,
                E const & _edges) const
    {
        if (points_.empty()) {
            _gnuplot << "print 'no point to process'\n;";
            return;
        }
        // pmin, pmax denotes bounding box
        point vmin = points_.front();
        point vmax = vmin;
        auto const pminmax = [&] (point const & p)
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
        std::for_each(std::cbegin(_vertices), std::cend(_vertices), [&] (auto const & v) { pminmax(v.c); });
        {
            _gnuplot << "set size square;\n"
                        "set key left;\n"
                        "unset colorbox;\n";
            _gnuplot << "set xrange [" << pmin.x << ':' << pmax.x << "];\n";
            _gnuplot << "set yrange [" << pmin.y << ':' << pmax.y << "];\n";
            _gnuplot << "set size ratio -1;\n";
        }
        {
            _gnuplot << "$sites << EOI\n";
            size_type i = 0;
            for (point const & point_ : points_) {
                _gnuplot << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
            }
            _gnuplot << "EOI\n";
        }
        if (draw_circles && !_vertices.empty()) {
            _gnuplot << "$circles << EOI\n";
            size_type i = 0;
            for (auto const & vertex_ : _vertices) {
                _gnuplot << vertex_.c.x << ' ' << vertex_.c.y << ' ' << vertex_.R << ' ' << i++ << '\n';
            }
            _gnuplot << "EOI\n";
        }
        auto const trunc_edge = [&] (point const & l, point const & r, point const & p) -> point
        {
            value_type const dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
            value_type const dy = l.x - r.x;
            auto const px = [&] (value_type const & y) -> point { return {(p.x + (y - p.y) * dx / dy), y}; };
            auto const py = [&] (value_type const & x) -> point
            {
                value_type const y = p.y + (x - p.x) * dy / dx;
                if (+eps < dy) {
                    if (vmax.y < y) {
                        return px(vmax.y);
                    }
                } else if (dy < -eps) {
                    if (y < vmin.y) {
                        return px(vmin.y);
                    }
                }
                return {x, y};
            };
            if (+eps < dx) {
                return py(vmax.x);
            } else if (dx < -eps) {
                return py(vmin.x);
            } else {
                if (+eps < dy) {
                    return {p.x, vmax.y};
                } else if (dy < -eps) {
                    return {p.x, vmin.y};
                } else {
                    assert(false);
                    return {p.x, p.y};
                }
            }
        };
        if (!_edges.empty()) {
            auto const pout = [&] (auto const & p)
            {
                _gnuplot << p.x << ' ' << p.y << '\n';
            };
            _gnuplot << "$edges << EOI\n";
            auto const nv = std::end(_vertices);
            for (auto const & edge_ : _edges) {
                bool const beg = (edge_.b != nv);
                bool const end = (edge_.e != nv);
                point const & l = *edge_.l;
                point const & r = *edge_.r;
                if (beg != end) {
                    assert(beg); // disable this, if flip at the mid of the sweepline<>::trunc_edge does not exist
                    point const & p = (beg ? edge_.b : edge_.e)->c;
                    if (!(p.x < vmin.x) && !(vmax.x < p.x) && !(p.y < vmin.y) && !(vmax.y < p.y)) {
                        pout(p);
                        pout(trunc_edge((beg ? l : r), (end ? l : r), p));
                    }
                } else if (beg) {
                    assert(!(edge_.e->c < edge_.b->c)); // disable this, if flip at the end of the sweepline<>::trunc_edge does not exist
                    pout(edge_.b->c);
                    pout(edge_.e->c);
                } else {
                    point const p{(l.x + r.x) / value_type(2), (l.y + r.y) / value_type(2)};
                    pout(trunc_edge(l, r, p));
                    pout(trunc_edge(r, l, p));
                }
                _gnuplot << "\n"; // separate lines
            }
            _gnuplot << "EOI\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '$sites' with points title 'sites # " << points_.size() << "'";
        if (draw_indices) {
            _gnuplot << ", '$sites' with labels offset character 0, character 1 notitle";
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

    void
    output(std::ostream & _out) const
    {
        return output(_out, sweepline_.vertices_, sweepline_.edges_);
    }

public :

    friend
    std::ostream &
    operator << (std::ostream & _gnuplot, voronoi const & _voronoi)
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
std::unique_ptr< char, decltype(std::free) & > demangled_name{nullptr, std::free};
std::size_t length = 0;
#pragma GCC diagnostic pop

}

inline
std::string
get_demangled_name(char const * const symbol) noexcept
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

    bool operator < (point const & p) const
    {
        return std::tie(x, y) < std::tie(p.x, p.y);
    }

};

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
        // Rectangle grid, diagonal grid or points uniformely distributed into a circle or square:
        {
            using seed_type = typename voronoi_type::seed_type;
#  if 0
            seed_type const seed = 2911579113;
#  else
            std::random_device D;
            auto const seed = static_cast< seed_type >(D());
#  endif
            voronoi_.seed(seed);
            //gnuplot_ << "set title 'seed = 0x" << std::hex << std::nouppercase << seed << ", N = " <<  std::dec << N << "'\n";
        }
#  if 0
        voronoi_.rectangle_grid(in_, 10); voronoi_.draw_circles = true;
#  elif 0
        voronoi_.diagonal_grid(in_, 20); voronoi_.draw_circles = true;
#  elif 0
        voronoi_.hexagonal_grid(in_, 20); voronoi_.eps = value_type(0.0001); //voronoi_.draw_circles = true;
#  elif 0
        voronoi_.triangular_grid(in_, 21); voronoi_.eps = value_type(0.0001); //voronoi_.draw_circles = true;
#  elif 0
        voronoi_.square(in_, value_type(10000), 100000);
#  else
        voronoi_.ball(in_, value_type(10000), 100000); // voronoi_.draw_circles = true; // voronoi_.draw_indices = true;
#  endif
# endif
        //log_ << in_.str() << '\n';
#endif
        {
            auto const start = steady_clock::now();
            in_ >> voronoi_;
            log_ << "input time = "
                 << duration_cast< microseconds >(steady_clock::now() - start).count()
                 << "us\n";
        }
        //voronoi_.swap_xy();
        //voronoi_.shift_xy(value_type(10000), value_type(10000));
    }
    { // run
        log_ << "start\n";
        auto const start = steady_clock::now();
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
        sweepline_type const & sweepline_ = voronoi_.sweepline_;
        log_ << "vertices # " << sweepline_.vertices_.size() << '\n';
        log_ << "edges # " << sweepline_.edges_.size() << '\n';
        std::string command_line_ = "sweepline.plt";
        std::ofstream f(command_line_);
        std::ostream & gnuplot_ = f;//std::cout;
#if 0
        { // clone (O(|vertices| * |edges|))
            using vertices = std::vector< typename sweepline_type::vertex >;
            using pvertex = typename vertices::const_iterator;
            using site = typename voronoi_type::site;
            vertices const vertices_{std::cbegin(sweepline_.vertices_), std::cend(sweepline_.vertices_)};
            pvertex const nv = std::cend(vertices_);
            auto const vclone = [&] (typename sweepline_type::pvertex const v) -> pvertex
            {
                return std::prev(nv, std::distance(v, sweepline_.nv));
            };
            struct edge { site l, r; pvertex b, e; };
            auto const eclone = [&] (typename sweepline_type::edge const & _edge) -> edge
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
#if 1
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
