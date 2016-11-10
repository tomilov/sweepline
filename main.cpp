#include "sweepline.hpp"

#include <utility>
#include <limits>
#include <chrono>
#include <iterator>
#include <random>
#include <chrono>
#include <tuple>
#include <vector>
#include <istream>
#include <ostream>

#include <cassert>
#include <cmath>

//#define SWEEPLINE_DRAW_CIRCLES
//#define SWEEPLINE_DRAW_INDICES

template< typename point_type, typename value_type >
struct voronoi
{

    using size_type = std::size_t;

    std::ostream & log_;

    // bounding box
    value_type const bbox = value_type(10000);

    value_type const eps = value_type(0.0001);

    value_type const delta = value_type(0.001);

    value_type const zero = value_type(0);
    value_type const one = value_type(1);

    voronoi(std::ostream & _log)
        : log_(_log)
    {
        assert(!(delta < eps));
        log_ << "eps = " << eps << '\n';
        log_ << "delta = " << delta << '\n';
    }

    using seed_type = typename std::mt19937::result_type;

    std::mt19937 rng;
    std::normal_distribution< value_type > normal_;
    std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};

    void
    seed(seed_type const seed)
    {
        log_ << "seed = " << seed << '\n';
        rng.seed(seed);
    }

    void
    uniform_circle(std::ostream & _out, size_type const N)
    {
        std::set< point_type, point_less > points_{point_less{delta}};
        _out << N << '\n';
        points_.clear();
        value_type const twosqreps = eps * (eps + eps);
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
            size_type m = 0;
            do {
                point_type p{normal_(rng), normal_(rng)};
                value_type norm = p.x * p.x + p.y * p.y;
                if (twosqreps < norm) {
                    using std::sqrt;
                    norm = bbox * sqrt(zero_to_one_(rng) / std::move(norm));
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

    struct ipoint { size_type x, y; };

    template< std::size_t nsqr >
    void
    quadrant(std::ostream & _out, size_type const max, ipoint const (& q)[nsqr])
    {
        if (0 == max) {
            _out << nsqr * 2 * 4 << '\n';
        } else {
            _out << (nsqr * 2 * 4) + 4 << '\n';
            _out << "0 " << max << '\n';
            _out << "0 -" << max << '\n';
            _out << max << " 0\n";
            _out << '-' << max << " 0\n";
        }
        auto const print = [&] (bool const direct, bool const signx, bool const signy)
        {
            for (ipoint const & p : q) {
                assert(p.x != 0);
                assert(p.y != 0);
                if (signx) {
                    _out << '-';
                }
                _out << (direct ? p.x : p.y) << ' ';
                if (signy) {
                    _out << '-';
                }
                _out << (direct ? p.y : p.x) << '\n';
            }
        };
        print(true,  true,  true);
        print(true,  false, true);
        print(true,  true,  false);
        print(true,  false, false);
        print(false, true,  true);
        print(false, false, true);
        print(false, true,  false);
        print(false, false, false);
    }

    using points = std::vector< point_type >;
    points sites_;

    using point_iterator = typename points::const_iterator;
    using sweepline_type = sweepline< point_iterator, point_type, value_type >;
    using point_less = typename sweepline_type::point_less;

    void input(std::istream & _in)
    {
        size_type N = 0;
        if (!(_in >> N)) {
            assert(false);
        }
        sites_.reserve(N);
        for (size_type n = 0; n < N; ++n) {
            point_type & point_ = sites_.emplace_back();
            if (!(_in >> point_.x)) {
                assert(false);
            }
            if (!(_in >> point_.y)) {
                assert(false);
            }
        }
        std::sort(std::begin(sites_), std::end(sites_), point_less{eps});
    }

    friend
    std::istream &
    operator >> (std::istream & _in, voronoi & _voronoi)
    {
        _voronoi.input(_in);
        return _in;
    }

    sweepline_type sweepline_{eps};

    void operator () ()
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        steady_clock::time_point const start = steady_clock::now();
        log_ << "begin sweepline\n";
        sweepline_(std::cbegin(sites_), std::cend(sites_));
        log_ << "sweepline time = "
             << duration_cast< microseconds >(steady_clock::now() - start).count()
             << "us\n";
    }

    value_type const vbox = value_type(1.5) * bbox;

    point_type
    trunc_edge(point_type const & l, point_type const & r, point_type const & p) const
    {
        value_type const dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
        value_type const dy = l.x - r.x;
        auto const px = [&] (value_type const & y) -> point_type { return {(p.x + (y - p.y) * dx / dy), y}; };
        auto const py = [&] (value_type const & x) -> point_type
        {
            value_type const y = p.y + (x - p.x) * dy / dx;
            if (+eps < dy) {
                if (+vbox < y) {
                    return px(+vbox);
                }
            } else if (dy < -eps) {
                if (y < -vbox) {
                    return px(-vbox);
                }
            }
            return {x, y};
        };
        if (+eps < dx) {
            return py(+vbox);
        } else if (dx < -eps) {
            return py(-vbox);
        } else {
            if (+eps < dy) {
                return {p.x, +vbox};
            } else if (dy < -eps) {
                return {p.x, -vbox};
            } else {
                assert(false);
                return {p.x, p.y};
            }
        }
    }

    void output(std::ostream & _gnuplot) const
    {
        if (sites_.empty()) {
            _gnuplot << "print 'no point to process'\n;";
            return;
        }
        {
            _gnuplot << "set size square;\n"
                        //"set size ratio -1;\n"
                        "set key left;\n"
                        "unset colorbox;\n";
            _gnuplot << "set xrange [" << -vbox << ':' << vbox << "];\n";
            _gnuplot << "set yrange [" << -vbox << ':' << vbox << "];\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '-' with points notitle";
#ifdef SWEEPLINE_DRAW_INDICES
        _gnuplot << ", '' with labels offset character 0, character 1 notitle";
#endif
#ifdef SWEEPLINE_DRAW_CIRCLES
        if (!sweepline_.vertices_.empty()) {
            _gnuplot << ", '' with circles notitle linecolor palette";
        }
#endif
        if (!sweepline_.edges_.empty()) {
            _gnuplot << ", '' with lines title 'edges (" << sweepline_.edges_.size() <<  ")'";
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
#ifdef SWEEPLINE_DRAW_INDICES
        {
            size_type i = 0;
            for (point_type const & point_ : sites_) {
                _gnuplot << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
            }
            _gnuplot << "e\n";
        }
#endif
#ifdef SWEEPLINE_DRAW_CIRCLES
        if (!sweepline_.vertices_.empty()) {
            size_type i = 0;
            for (auto const & vertex_ : sweepline_.vertices_) {
                _gnuplot << vertex_.c.x << ' ' << vertex_.c.y << ' ' << vertex_.R << ' ' << i++ << '\n';
            }
            _gnuplot << "e\n";
        }
#endif
        if (!sweepline_.edges_.empty()) {
            for (auto const & edge_ : sweepline_.edges_) {
                bool const beg = (edge_.b != sweepline_.nov);
                bool const end = (edge_.e != sweepline_.nov);
                point_type const & l = *edge_.l;
                point_type const & r = *edge_.r;
                if (beg != end) {
                    point_type const & p = (beg ? edge_.b : edge_.e)->c;
                    if (!(p.x < -vbox) && !(vbox < p.x) && !(p.y < -vbox) && !(vbox < p.y)) {
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

    friend
    std::ostream &
    operator << (std::ostream & _gnuplot, voronoi const & _voronoi)
    {
        _voronoi.output(_gnuplot);
        return _gnuplot;
    }

};

#include <iomanip>
#include <iostream>
#include <sstream>

#include <cstdlib>

#include <x86intrin.h>

using value_type = double;

struct point
{

    value_type x, y;

};

using point_type = point;

int main()
{
    using voronoi_type = voronoi< point_type, value_type >;
    voronoi_type voronoi_{std::clog};
    std::ostream & gnuplot_ = std::cout;
    {
#if 0
        std::istream & in_ = std::cin;
#elif 1
        std::stringstream in_;
        in_ >> std::scientific;
        in_.precision(std::numeric_limits< value_type >::digits10 + 2);
#if 1
        in_ << "7\n"
               "1 0\n"
               "2 0\n"
               "3 0\n"
               "4 0\n"
               "5 0\n"
               "6 0\n"
               "7 0\n";
#elif 0
        // Concentric:
#if 0
        in_ << "4\n"
               "-1 0\n"
               "0 -1\n"
               "0 1\n"
               "1 0\n";
#elif 0
        voronoi_.quadrant(in_, 5, {{3, 4}});
#elif 0
        voronoi_.quadrant(in_, 25, {{7, 24}, {15, 20}});
#elif 0
        voronoi_.quadrant(in_, 65, {{16, 63}, {25, 60}, {33, 56}, {39, 52}});
#elif 0
        voronoi_.quadrant(in_, 325, {{36, 323}, {80, 315}, {91, 312}, {125, 300}, {165, 280}, {195, 260}, {204, 253}});
#elif 0
        voronoi_.quadrant(in_, 1105, {{47,  1104}, {105, 1100}, {169, 1092}, {264, 1073}, {272, 1071}, {425, 1020}, {468, 1001},
                                      {520, 975 }, {561, 952 }, {576, 943 }, {663, 884 }, {700, 855 }, {744, 817 }});
#else
        voronoi_.quadrant(in_, 5525, {{235,  5520}, {525,  5500}, {612,  5491}, {845,  5460},
                                      {1036, 5427}, {1131, 5408}, {1320, 5365}, {1360, 5355}, {1547, 5304},
                                      {2044, 5133}, {2125, 5100}, {2163, 5084}, {2340, 5005}, {2600, 4875},
                                      {2805, 4760}, {2880, 4715}, {3124, 4557}, {3315, 4420}, {3468, 4301},
                                      {3500, 4275}, {3720, 4085}, {3861, 3952}});
#endif
#elif 0
        // Uniformely distributed into the circle
        constexpr std::size_t N = 100000;
        {
            using seed_type = typename voronoi_type::seed_type;
#if 0
            seed_type const seed = 2847645394;
#else
            std::random_device D;
            auto const seed = static_cast< seed_type >(D());
#endif
            voronoi_.seed(seed);
            gnuplot_ << "set title 'seed = 0x" << std::hex << seed << ", N = " <<  std::dec << N << "'\n";
        }
        voronoi_.uniform_circle(in_, N);
#endif
        std::clog << in_.str() << '\n';
#endif
        in_ >> voronoi_;
    }
    voronoi_();
    gnuplot_ << voronoi_ << std::endl;
    return EXIT_SUCCESS;
}
