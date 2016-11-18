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

template< typename point_type, typename value_type = decltype(std::declval< point_type >().x) >
struct voronoi
{

    using size_type = std::size_t;

    std::ostream & log_;

    bool draw_indices = false;
    bool draw_circles = false;

    // bounding box

    value_type eps = value_type(10) * std::numeric_limits< value_type >::epsilon();

    value_type delta = value_type(0.001);

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

    void
    ball(std::ostream & _out, value_type const radius, size_type const N)
    {
        std::set< point_type, point_less > points_{point_less{delta}};
        _out << N << '\n';
        value_type const twosqreps = eps * (eps + eps);
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
            size_type m = 0;
            do {
                point_type p{normal_(rng), normal_(rng)};
                value_type norm = p.x * p.x + p.y * p.y;
                if (twosqreps < norm) {
                    using std::sqrt;
                    norm = radius * sqrt(zero_to_one_(rng) / std::move(norm));
                    p.x *= norm;
                    p.y *= norm;
                } else {
                    p.x = p.y = zero;
                }
                if (points_.insert(std::move(p)).second) {
                    break;
                }
            } while (++m < M);
            if (m == M) {
                log_ << "the number (" << M << ") of attempts is exceeded\n";
                log_ << "only " << n << "points generated\n";
                break;
            }
        }
        for (point_type const & point_ : points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    void
    square(std::ostream & _out, value_type const bbox, size_type const N)
    {
        std::set< point_type, point_less > points_{point_less{delta}};
        _out << N << '\n';
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed square
            size_type m = 0;
            do {
                point_type p{zero_to_one_(rng), zero_to_one_(rng)};
                p.x += p.x;
                p.y += p.y;
                p.x -= one;
                p.y -= one;
                p.x *= bbox;
                p.y *= bbox;
                if (points_.insert(std::move(p)).second) {
                    break;
                }
            } while (++m < M);
            if (m == M) {
                log_ << "the number (" << M << ") of attempts is exceeded\n";
                log_ << "only " << n << "points generated\n";
                break;
            }
        }
        for (point_type const & point_ : points_) {
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
        size_type const N = 2 * (size * (size + 1) - (size % 2));
        _out << N << '\n';
        using std::sqrt;
        value_type const step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            value_type const xx = x * step;
            if ((x % 2) == 0) {
                value_type const yy = 3 * (x - 1);
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
        assert(i == N);
    }

    void
    triangular_grid(std::ostream & _out, size_type const size) const
    {
        size_type const N = 2 * (1 + size + size) * size;
        _out << N << '\n';
        using std::sqrt;
        value_type const step = sqrt(value_type(3));
        size_type i = 0;
        for (size_type x = 1; x <= size; ++x) {
            size_type const xx = 3 * x;
            {
                value_type const yy = step * (xx - 1 - (x % 2));
                _out << "0 " << yy << '\n';
                _out << "0 -" << yy << '\n';
            }
            i += 2;
            for (size_type y = 1; y <= size; ++y) {
                value_type const yy = step * (3 * y - 1 - ((y + x) % 2));
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

    using points = std::vector< point_type >;
    using site = typename points::const_iterator;

    using sweepline_type = sweepline< site, point_type, value_type >;

private :

    using point_less = typename sweepline_type::point_less;

    points sites_;

    void input(std::istream & _in)
    {
        size_type N = 0;
        if (!(_in >> N)) {
            assert(false);
        }
        sites_.reserve(N);
        for (size_type n = 0; n < N; ++n) {
            sites_.emplace_back();
            point_type & point_ = sites_.back();
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
        for (point_type & point_ : sites_) {
            using std::swap;
            swap(point_.x, point_.y);
        }
    }

    void
    shift_xy(value_type const & dx, value_type const & dy)
    {
        for (point_type & point_ : sites_) {
            point_.x += dx;
            point_.y += dy;
        }
    }

    sweepline_type sweepline_{eps};

    void operator () ()
    {
        assert((std::set< point_type, point_less >{std::cbegin(sites_), std::cend(sites_), point_less{eps}}.size() == sites_.size()));
        log_ << "N = " << sites_.size() << '\n';
#if 0
        using pproxy = std::vector< site >;
        pproxy pproxy_;
        pproxy_.reserve(sites_.size());
        auto const send = std::cend(sites_);
        for (auto p = std::cbegin(sites_); p != send; ++p) {
            pproxy_.push_back(p);
        }
        std::sort(std::begin(pproxy_), std::end(pproxy_), point_less{eps});
        using ppoint = typename pproxy::const_iterator;
        struct point_proxy
                : std::iterator< std::forward_iterator_tag, point_type const >
        {

            ppoint p;

            point_proxy(ppoint pp) : p(pp) { ; }
            operator site const & () const { return *p; }

            point_proxy & operator ++ () { ++p; return *this; }
            point_proxy operator ++ (int) { auto pp = p; operator ++ (); return {pp}; }

            point_type const & operator * () const { return **p; }

            bool operator == (point_proxy const & rhs) const { return (p == rhs.p); }
            bool operator != (point_proxy const & rhs) const { return !operator == (rhs); }

        };
        sweepline_(point_proxy{std::cbegin(pproxy_)}, point_proxy{std::cend(pproxy_)});
#else
        std::sort(std::begin(sites_), std::end(sites_), point_less{eps});
        sweepline_(std::cbegin(sites_), std::cend(sites_));
#endif
    }

    value_type zoom = value_type(0.2);

    template< typename V, typename E >
    void output(std::ostream & _gnuplot,
                V const & _vertices,
                E const & _edges) const
    {
        if (sites_.empty()) {
            _gnuplot << "print 'no point to process'\n;";
            return;
        }
        // pmin, pmax denotes bounding box
        point_type vmin = sites_.front();
        point_type vmax = vmin;
        auto const pminmax = [&] (point_type const & p)
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
        std::for_each(std::next(std::cbegin(sites_)), std::cend(sites_), pminmax);
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
        point_type pmin = vmin;
        point_type pmax = vmax;
        std::for_each(std::cbegin(_vertices), std::cend(_vertices), [&] (auto const & v) { pminmax(v.c); });
        {
            _gnuplot << "set size square;\n"
                        "set key left;\n"
                        "unset colorbox;\n";
            _gnuplot << "set xrange [" << pmin.x << ':' << pmax.x << "];\n";
            _gnuplot << "set yrange [" << pmin.y << ':' << pmax.y << "];\n";
            _gnuplot << "set size ratio -1;\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '-' with points title 'sites # " << sites_.size() << "'";
        if (draw_indices) {
            _gnuplot << ", '' with labels offset character 0, character 1 notitle";
        }
        if (draw_circles && !_vertices.empty()) {
            _gnuplot << ", '' with circles title 'vertices # " << _vertices.size() << "' linecolor palette";
        }
        if (!_edges.empty()) {
            _gnuplot << ", '' with lines title 'edges # " << _edges.size() <<  "'";
        }
        _gnuplot << ";\n";
        auto const pout = [&] (auto const & p)
        {
            _gnuplot << p.x << ' ' << p.y << '\n';
        };
        {
            for (point_type const & point_ : sites_) {
                pout(point_);
            }
            _gnuplot << "e\n";
        }
        if (draw_indices) {
            size_type i = 0;
            for (point_type const & point_ : sites_) {
                _gnuplot << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
            }
            _gnuplot << "e\n";
        }
        if (draw_circles && !_vertices.empty()) {
            size_type i = 0;
            for (auto const & vertex_ : _vertices) {
                _gnuplot << vertex_.c.x << ' ' << vertex_.c.y << ' ' << vertex_.R << ' ' << i++ << '\n';
            }
            _gnuplot << "e\n";
        }
        auto const trunc_edge = [&] (point_type const & l, point_type const & r, point_type const & p) -> point_type
        {
            value_type const dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
            value_type const dy = l.x - r.x;
            auto const px = [&] (value_type const & y) -> point_type { return {(p.x + (y - p.y) * dx / dy), y}; };
            auto const py = [&] (value_type const & x) -> point_type
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
            auto const nov = std::end(_vertices);
            for (auto const & edge_ : _edges) {
                bool const beg = (edge_.b != nov);
                bool const end = (edge_.e != nov);
                point_type const & l = *edge_.l;
                point_type const & r = *edge_.r;
                if (beg != end) {
                    point_type const & p = (beg ? edge_.b : edge_.e)->c;
                    if (!(p.x < vmin.x) && !(vmax.x < p.x) && !(p.y < vmin.y) && !(vmax.y < p.y)) {
                        pout(p);
                        pout(trunc_edge((beg ? l : r), (end ? l : r), p));
                        _gnuplot << "\n";
                    }
                } else if (beg && end) {
                    pout(edge_.b->c);
                    pout(edge_.e->c);
                    _gnuplot << "\n";
                } else {
                    point_type const p{(l.x + r.x) / value_type(2), (l.y + r.y) / value_type(2)};
                    pout(trunc_edge(l, r, p));
                    pout(trunc_edge(r, l, p));
                }
                _gnuplot << "\n"; // separate lines (sic! it is incredible, but output for chained lines works fine!!!)
            }
            _gnuplot << "e\n";
        }
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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
std::mutex m;
std::unique_ptr< char, decltype(std::free) & > demangled_name{nullptr, std::free};
std::size_t length = 0;
#pragma clang diagnostic pop

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

};

int main()
{
    using voronoi_type = voronoi< point >;
    std::ostream & log_ = std::clog;
    voronoi_type voronoi_{log_};
    { // input
#if 0
        std::istream & in_ = std::cin;
#else
        std::stringstream in_;
        in_ >> std::hexfloat;
        in_.precision(std::numeric_limits< value_type >::digits10 + 2);
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
# elif 1
        // Rectangle grid, diagonal grid or points uniformely distributed into a circle or square:
        {
            using seed_type = typename voronoi_type::seed_type;
#  if 0
            seed_type const seed = 855215359;
#  else
            std::random_device D;
            auto const seed = static_cast< seed_type >(D());
#  endif
            voronoi_.seed(seed);
            //gnuplot_ << "set title 'seed = 0x" << std::hex << std::nouppercase << seed << ", N = " <<  std::dec << N << "'\n";
        }
        //voronoi_.rectangle_grid(in_, 10); voronoi_.draw_circles = true;
        //voronoi_.diagonal_grid(in_, 20); voronoi_.draw_circles = true;
        //voronoi_.hexagonal_grid(in_, 20); //voronoi_.draw_circles = true;
        //voronoi_.triangular_grid(in_, 200); voronoi_.eps = value_type(0.0001); //voronoi_.draw_circles = true;
        voronoi_.ball(in_, value_type(10000), 15); voronoi_.draw_circles = true; // voronoi_.draw_indices = true;
        //voronoi_.square(in_, value_type(10000), 100000);
# endif
        //log_ << in_.str() << '\n';
#endif
        in_ >> voronoi_;
        //voronoi_.swap_xy();
        //voronoi_.shift_xy(value_type(10000), value_type(10000));
    }
    { // run
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
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
        log_ << "begin sweepline\n";
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
        std::ostream & gnuplot_ = std::cout;
#if 0
        { // clone
            using vertices = std::vector< typename sweepline_type::vertex >;
            using pvertex = typename vertices::const_iterator;
            using site = typename voronoi_type::site;
            vertices const vertices_{std::cbegin(sweepline_.vertices_), std::cend(sweepline_.vertices_)};
            pvertex const nov = std::cend(vertices_);
            auto const vclone = [&] (typename sweepline_type::pvertex const v) -> pvertex
            {
                return std::prev(nov, std::distance(v, sweepline_.nov));
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
    }
    return EXIT_SUCCESS;
}
