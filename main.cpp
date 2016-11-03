#include "sweepline.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <vector>
#include <chrono>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <x86intrin.h>

template< typename point_type, typename point_less, typename value_type >
struct voronoi
{

    using size_type = std::size_t;

    std::ostream & log_;

    value_type const eps = [] { return value_type(1E-9); }();
    value_type const zero = value_type(0);
    value_type const one = value_type(1);

    // bounding box
    value_type const bbox = value_type(10);
    value_type const delta = eps * value_type(10);

    std::ostream & gnuplot_ = std::cout;

    void
    generate(std::ostream & _out, size_type const N = 100000)
    {
        using seed_type = typename std::mt19937::result_type;
#if 1
        // ss == 953, 934 seed = 0x13d69d450e99 N == 1000
        seed_type const seed_ = 15959779189989;
#elif 0
        std::random_device D;
        auto const seed_ = static_cast< seed_type >(D());
#elif 1
        //auto const seed_ = static_cast< seed_type >(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        auto const seed_ = static_cast< seed_type >(__rdtsc());
#endif
        log_ << seed_ << '\n';
        gnuplot_ << "set title 'seed = 0x" << std::hex << seed_ << ", N = " <<  std::dec << N << "'\n";
        std::mt19937 g{seed_};
        std::normal_distribution< value_type > normal_;
        std::set< point_type, point_less > points_{point_less{delta}};
        std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};
        _out << N << "\n";
        points_.clear();
        value_type const twosqreps = eps * (eps + eps);
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
            for (;;) {
                point_type p{normal_(g), normal_(g)};
                value_type norm = p.x * p.x + p.y * p.y;
                if (twosqreps < norm) {
                    using std::sqrt;
                    norm = bbox * sqrt(zero_to_one_(g) / std::move(norm));
                    p.x *= norm;
                    p.y *= norm;
                } else {
                    p.x = p.y = zero;
                }
                if (points_.insert(std::move(p)).second) {
                    break;
                }
            }
        }
        for (point_type const & point_ : points_) {
            _out << point_.x << ' ' << point_.y << '\n';
        }
    }

    voronoi(std::ostream & _log = std::clog)
        : log_(_log)
    {
        assert(!(delta < eps));
    }

    using points = std::vector< point_type >;
    points sites_;

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

    using point_iterator = typename points::const_iterator;
    using sweepline_type = sweepline< point_iterator, point_less, point_type, value_type >;
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

    value_type const vbox = value_type(2) * bbox;

    point_type
    trunc_edge(point_type const & l, point_type const & r, point_type const & p) const
    {
        value_type const dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
        value_type const dy = l.x - r.x;
        auto const pp = [&] (value_type const & y) -> point_type { return {(p.x + (y - p.y) * dx / dy), y}; };
        auto const yy = [&] (value_type const & x) { return p.y + (x - p.x) * dy / dx; };
        if (eps < dx) {
            value_type const y = yy(vbox);
            if (eps < dy) {
                if (vbox < y) {
                    return pp(vbox);
                }
            } else if (dy < -eps) {
                if (y < -vbox) {
                    return pp(-vbox);
                }
            }
            return {vbox, y};
        } else if (dx < eps) {
            value_type const y = yy(-vbox);
            if (eps < dy) {
                if (vbox < y) {
                    return pp(vbox);
                }
            } else if (dy < -eps) {
                if (y < -vbox) {
                    return pp(-vbox);
                }
            }
            return {-vbox, y};
        } else {
            if (eps < dy) {
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
                        "set key left;\n";
            _gnuplot << "set xrange [" << -vbox << ':' << vbox << "];\n";
            _gnuplot << "set yrange [" << -vbox << ':' << vbox << "];\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '-' with points notitle"
                    ", '' with labels offset character 0, character 1 notitle";
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
        {
            size_type i = 0;
            for (point_type const & point_ : sites_) {
                _gnuplot << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
            }
            _gnuplot << "e\n";
        }
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

using value_type = double;

struct point
{

    value_type x, y;

    friend
    std::ostream &
    operator << (std::ostream & _out, point const & p)
    {
        return _out << '{' << p.x << ", " << p.y << '}';
    }

};

using point_type = point;

struct point_less
{

    value_type const & eps_;

    bool operator () (point_type const & _lhs, point_type const & _rhs) const
    {
        value_type const & x = _lhs.x + eps_;
        value_type const & y = _lhs.y + eps_;
        return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
    }

};

int
main()
{
    voronoi< point_type, point_less, value_type > voronoi_{std::clog};
    {
#if 0
        std::istream & in_ = std::cin;
#elif 1
        std::stringstream in_;
        in_ >> std::scientific;
        in_.precision(std::numeric_limits< value_type >::digits10 + 2);
        voronoi_.generate(in_, 1000);
        std::clog << in_.str() << '\n';
#elif 0
        std::stringstream in_;
        in_ << "3\n"
               "-1 0\n"
               "0 -1\n"
               "0 1\n";
#endif
        in_ >> voronoi_;
    }
    voronoi_();
    std::cout << voronoi_ << std::endl;
    return EXIT_SUCCESS;
}
