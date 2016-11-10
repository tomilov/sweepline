#include "sweepline.hpp"

#include <utility>
#include <limits>
#include <chrono>
#include <iterator>
#include <algorithm>
#include <random>
#include <tuple>
#include <set>
#include <vector>
#include <iomanip>
#include <istream>
#include <ostream>

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

    value_type eps = value_type(1E-14);

    value_type delta = value_type(0.001);

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

private :

    std::mt19937 rng;
    std::normal_distribution< value_type > normal_;
    std::uniform_real_distribution< value_type > zero_to_one_{zero, std::nextafter(one, one + one)};

public :

    void
    seed(seed_type const seed)
    {
        log_ << "seed = " << seed << '\n';
        rng.seed(seed);
    }

    void
    uniform_circle(std::ostream & _out, value_type const bbox, size_type const N)
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

    void
    uniform_square(std::ostream & _out, value_type const bbox, size_type const N)
    {
        std::set< point_type, point_less > points_{point_less{delta}};
        _out << N << '\n';
        constexpr size_type M = 1000; // number of attempts
        for (size_type n = 0; n < N; ++n) { // points that are uniformely distributed inside of closed ball
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
    rectangle_mesh(std::ostream & _out, size_type const bbox) const
    {
        _out << (1 + 4 * bbox * (bbox + 1)) << '\n';
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

    using sweepline_type = sweepline< site >;

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

public :

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
        sweepline_(std::cbegin(sites_), std::cend(sites_));
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
        point_type pmin = sites_.front();
        point_type pmax = pmin;
        auto const minmax = [&] (point_type const & p)
        {
            if (p.x < pmin.x) {
                pmin.x = p.x;
            } else if (pmax.x < p.x) {
                pmax.x = p.x;
            }
            if (p.y < pmin.y) {
                pmin.y = p.y;
            } else if (pmax.y < p.y) {
                pmax.y = p.y;
            }
        };
        std::for_each(std::next(std::cbegin(sites_)), std::cend(sites_), minmax);
        assert(value_type(-0.5) < zoom);
        if (pmin.x + delta < pmax.x) {
            value_type dx = pmax.x - pmin.x;
            dx *= zoom;
            pmin.x -= dx;
            pmax.x += dx;
        } else {
            pmin.x -= delta;
            pmax.x += delta;
        }
        if (pmin.y + delta < pmax.y) {
            value_type dy = pmax.y - pmin.y;
            dy *= zoom;
            pmin.y -= dy;
            pmax.y += dy;
        } else {
            pmin.y -= delta;
            pmax.y += delta;
        }
        {
            _gnuplot << "set size square;\n"
                        "set key left;\n"
                        "unset colorbox;\n";
            if (draw_circles) {
                _gnuplot << "set size ratio -1;\n";
            }
            _gnuplot << "set xrange [" << pmin.x << ':' << pmax.x << "];\n";
            _gnuplot << "set yrange [" << pmin.y << ':' << pmax.y << "];\n";
        }
        _gnuplot << "plot";
        _gnuplot << " '-' with points notitle";
        if (draw_indices) {
            _gnuplot << ", '' with labels offset character 0, character 1 notitle";
        }
        if (draw_circles && !_vertices.empty()) {
            _gnuplot << ", '' with circles notitle linecolor palette";
        }
        if (!_edges.empty()) {
            _gnuplot << ", '' with lines title 'edges (" << _edges.size() <<  ")'";
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
                    if (pmax.y < y) {
                        return px(pmax.y);
                    }
                } else if (dy < -eps) {
                    if (y < pmin.y) {
                        return px(pmin.y);
                    }
                }
                return {x, y};
            };
            if (+eps < dx) {
                return py(pmax.x);
            } else if (dx < -eps) {
                return py(pmin.x);
            } else {
                if (+eps < dy) {
                    return {p.x, pmax.y};
                } else if (dy < -eps) {
                    return {p.x, pmin.y};
                } else {
                    assert(false);
                    return {p.x, p.y};
                }
            }
        };
        if (!_edges.empty()) {
            for (auto const & edge_ : _edges) {
                bool const beg = (edge_.b != std::end(_vertices));
                bool const end = (edge_.e != std::end(_vertices));
                point_type const & l = *edge_.l;
                point_type const & r = *edge_.r;
                if (beg != end) {
                    point_type const & p = (beg ? edge_.b : edge_.e)->c;
                    if (!(p.x < pmin.x) && !(pmax.x < p.x) && !(p.y < pmin.y) && !(pmax.y < p.y)) {
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

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <cstdlib>

using value_type = double;

struct point
{

    value_type x, y;

};

int main()
{
    using voronoi_type = voronoi< point >;
    std::ostream & log_ = std::clog;
    voronoi_type voronoi_{log_};
    // input:
    std::ostream & gnuplot_ = std::cout;
    {
#if 0
        std::istream & in_ = std::cin;
#elif 1
        std::stringstream in_;
        in_ >> std::scientific;
        in_.precision(std::numeric_limits< value_type >::digits10 + 2);
#if 0
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
#elif 1
        // Uniformely distributed into the circle or square
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
            gnuplot_ << "set title 'seed = 0x" << std::hex << std::nouppercase << seed << ", N = " <<  std::dec << N << "'\n";
        }
        //voronoi_.uniform_circle(in_, value_type(10000), N);
        //voronoi_.uniform_square(in_, value_type(10000), N);
        voronoi_.rectangle_mesh(in_, 10);
#endif
        //log_ << in_.str() << '\n';
#endif
        in_ >> voronoi_;
    }
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::steady_clock;
        auto const start = steady_clock::now();
        voronoi_();
        log_ << "begin sweepline\n";
        log_ << "sweepline time = "
             << duration_cast< microseconds >(steady_clock::now() - start).count()
             << "us\n";
    }
    // output:
    // setup:
    voronoi_.draw_circles = false; // (sweepline_.vertices_.size() < 300);
    voronoi_.draw_indices = false;
    auto const & sweepline_ = voronoi_.sweepline_;
    log_ << "vertices # " << sweepline_.vertices_.size() << '\n';
    log_ << "edges # " << sweepline_.edges_.size() << '\n';
#if 1
    gnuplot_ << voronoi_ << std::endl;
#else
    { // clone
        using sweepline_type = typename voronoi_type::sweepline_type;
        // source
        // destination
        using vertices = std::vector< typename sweepline_type::vertex >;
        using pvertex = typename vertices::const_iterator;
        using site = typename voronoi_type::site;
        // transformation function:
        vertices const vertices_{std::cbegin(sweepline_.vertices_), std::cend(sweepline_.vertices_)};
        pvertex const nov = std::cend(vertices_);
        // cloning
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
#endif
    return EXIT_SUCCESS;
}
