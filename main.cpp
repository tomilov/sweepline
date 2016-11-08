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

template< typename point_type, typename point_less, typename value_type >
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
                        "set size ratio -1;\n"
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

};

using point_type = point;

struct point_less
{

    value_type const & eps_;

    bool operator () (point_type const & _lhs, point_type const & _rhs) const
    {
        if (_lhs.x + eps_ < _rhs.x) {
            return true;
        } else if (_rhs.x + eps_ < _lhs.x) {
            return false;
        } else {
            if (_lhs.y + eps_ < _rhs.y) {
                return true;
            } else {
                return false;
            }
        }
    }

};

#include <iomanip>
#include <iostream>
#include <sstream>

#include <cstdlib>

#include <x86intrin.h>

int
main()
{
    constexpr std::size_t N = 100000;

    using voronoi_type = voronoi< point_type, point_less, value_type >;
    voronoi_type voronoi_{std::clog};
    std::ostream & gnuplot_ = std::cout;
    {
#if 0
        std::istream & in_ = std::cin;
#else
        std::stringstream in_;
        in_ >> std::scientific;
        in_.precision(std::numeric_limits< value_type >::digits10 + 2);
        {
            using seed_type = typename voronoi_type::seed_type;
#if 1
            seed_type const seed = 2847645394;
#else
            std::random_device D;
            auto const seed = static_cast< seed_type >(D());
#endif
            voronoi_.seed(seed);
            gnuplot_ << "set title 'seed = 0x" << std::hex << seed << ", N = " <<  std::dec << N << "'\n";
        }
        voronoi_.uniform_circle(in_, N);
        //std::clog << in_.str() << '\n';
#endif
        in_ >> voronoi_;
    }
    voronoi_();
    gnuplot_ << voronoi_ << std::endl;
    return EXIT_SUCCESS;
}
