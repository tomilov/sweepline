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

using size_type = std::size_t;

using value_type = double;

struct point
{

    value_type x, y;

    bool operator < (point const & _rhs) const
    {
        return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
    }

};

struct point_less
{

    value_type const & eps_;

    bool operator () (point const & _lhs, point const & _rhs) const
    {
        value_type const & x = _lhs.x + eps_;
        value_type const & y = _lhs.y + eps_;
        return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
    }

};

namespace
{

// bounding box
value_type const bbox = value_type(100);
value_type const delta = value_type(1E-12);

value_type const eps = std::numeric_limits< value_type >::epsilon();
value_type const zero = value_type(0);
value_type const one = value_type(1);

point_less const point_less_{delta};

std::ostream & gnuplot_ = std::cout;
std::ostream & log_ = std::clog;

void
generate(std::ostream & _out, size_type const N = 100000)
{
    using seed_type = typename std::mt19937::result_type;
#if 0
    seed_type const seed_ = 23902348445254;//19614518643971;//8864935383105;
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
    std::set< point, point_less > points_{point_less_};
    std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};
    _out << N << "\n";
    points_.clear();
    value_type const twosqreps = eps * (eps + eps);
    for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
        for (;;) {
            point p{normal_(g), normal_(g)};
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
    for (point const & point_ : points_) {
        _out << point_.x << ' ' << point_.y << '\n';
    }
}

} // namespace

int
main()
{
    assert(!(delta < eps));
#if 0
    std::istream & in_ = std::cin;
#elif 1
    std::stringstream in_;
    generate(in_, 100000);
#elif 0
    std::stringstream in_;
    in_ << "3\n"
           "0 0\n"
           "1 -0.1\n"
           "3 -1\n";
#elif 0
    std::stringstream in_;
    in_ << "2\n"
           "1 0\n"
           "0 1\n";
#endif
    size_type N{};
    if (!(in_ >> N)) {
        assert(false);
    }
    assert(0 < N);
    using point_type = point;
    using points = std::vector< point_type >;
    points points_;
    points_.reserve(N);
    for (size_type n = 0; n < N; ++n) {
        point_type & point_ = points_.emplace_back();
        if (!(in_ >> point_.x)) {
            assert(false);
        }
        if (!(in_ >> point_.y)) {
            assert(false);
        }
    }
    using point_iterator = typename points::const_iterator;
    using sweepline_type = sweepline< point_iterator, point_type const, value_type >;
    sweepline_type sweepline_{eps};
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        steady_clock::time_point const start = steady_clock::now();
        log_ << "begin sweepline\n";
        sweepline_(std::cbegin(points_), std::cend(points_));
        log_ << "sweepline time = "
             << duration_cast< microseconds >(steady_clock::now() - start).count()
             << "us" << std::endl;
    }
    {
        value_type const vbox = value_type(2) * bbox;
        {
            gnuplot_ << "set size square;\n"
                        "set key left;\n";
            gnuplot_ << "set xrange [" << -vbox << ':' << vbox << "];\n";
            gnuplot_ << "set yrange [" << -vbox << ':' << vbox << "];\n";
        }
        {
            gnuplot_ << "plot";
            gnuplot_ << " '-' with points notitle"
                        ", '' with labels offset character 0, character 1 notitle";
            if (!sweepline_.edges_.empty()) {
                gnuplot_ << ", '' with lines title 'edges (" << sweepline_.edges_.size() <<  ")'";
            }
            gnuplot_ << ";\n";
            auto const pout = [&] (value_type const & x, value_type const & y)
            {
                gnuplot_ << x << ' ' << y << '\n';
            };
            {
                for (auto const & point_ : points_) {
                    pout(point_.x, point_.y);
                }
                gnuplot_ << "e\n";
            }
            {
                size_type i = 0;
                for (auto const & point_ : points_) {
                    gnuplot_ << point_.x << ' ' << point_.y << ' ' << i++ << '\n';
                }
                gnuplot_ << "e\n";
            }
            if (!sweepline_.edges_.empty()) {
                for (auto const & edge_ : sweepline_.edges_) {
                    auto const & l = *edge_.l;
                    auto const & r = *edge_.r;
                    value_type const dx = r.y - l.y; // +pi/2 rotation (dy, -dx)
                    value_type const dy = l.x - r.x;
                    auto const pend = [&] (auto const & p) -> point_type
                    {
                        if (eps < dx) {
                            value_type const yy = p.y + (vbox - p.x) * dy / dx;
                            if (eps < dy) {
                                if (vbox < yy) {
                                    value_type const xx = p.x + (vbox - p.y) * dx / dy;
                                    return {xx, vbox};
                                }
                            } else if (dy < -eps) {
                                if (yy < -vbox) {
                                    value_type const xx = p.x - (vbox + p.y) * dx / dy;
                                    return {xx, -vbox};
                                }
                            }
                            return {vbox, yy};
                        } else if (dx < eps) {
                            value_type const yy = p.y - (vbox + p.x) * dy / dx;
                            if (eps < dy) {
                                if (vbox < yy) {
                                    value_type const xx = p.x + (vbox - p.y) * dx / dy;
                                    return {xx, vbox};
                                }
                            } else if (dy < -eps) {
                                if (yy < -vbox) {
                                    value_type const xx = p.x - (vbox + p.y) * dx / dy;
                                    return {xx, -vbox};
                                }
                            }
                            return {-vbox, yy};
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
                    };
                    bool const beg = (edge_.b != sweepline_.nov);
                    bool const end = (edge_.e != sweepline_.nov);
                    if (beg && !end) {
                        auto const & p = *edge_.b;
                        if (!(p.x < -vbox) && !(vbox < p.x) && !(p.y < -vbox) && !(vbox < p.y)) {
                            pout(p.x, p.y);
                            auto const & e = pend(p);
                            pout(e.x, e.y);
                            gnuplot_ << "\n";
                        }
                    } else if (beg && end) {
                        auto const & b = *edge_.b;
                        pout(b.x, b.y);
                        auto const & e = *edge_.e;
                        pout(e.x, e.y);
                        gnuplot_ << "\n";
                    } else {
                        assert(!beg && !end);
                        //value_type const x = (l.x + r.x) / value_type(2);
                        //value_type const y = (l.y + r.y) / value_type(2);
                        assert(false && "need to implement");
                    }
                }
                gnuplot_ << "e\n";
            }
        }
    }
    return EXIT_SUCCESS;
}
