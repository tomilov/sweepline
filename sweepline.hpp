#pragma once

#include <algorithm>
#include <deque>
#include <functional>
#include <numeric>
#include <iterator>
#include <limits>
#include <list>
#include <set>
#include <map>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <experimental/optional>
#include <stdexcept>

#include <cassert>
#include <cmath>

template< typename point_iterator,
          typename point_type, // lexicagraphically comparable (x then y)
          typename value_type >
struct sweepline
{

    static_assert(std::is_same_v< decltype(std::declval< point_type >().x), decltype(std::declval< point_type >().y) >, "point_type should have x and y data members of the same type");

    using size_type = std::size_t;

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct site
    {

        site(point_iterator _p)
            : p(std::move(_p))
        { ; }

        operator point_iterator () const
        {
            return p;
        }

        point_iterator
        operator -> () const
        {
            return p;
        }

        point_type &
        operator * () const
        {
            return *p;
        }

        bool operator < (site const & _rhs) const
        {
            return *p < *_rhs;
        }

    private :

        point_iterator p;

    };

    struct vertex // denote circumscribed circle
    {

        point_type p; // circumcenter
        value_type R; // radius

        value_type x() const { return p.x + R; }
        value_type const & y() const { return p.y; }

    };

    struct vertex_less
    {

        value_type const & eps_;

        bool operator () (vertex const & _lhs, vertex const & _rhs) const // lexicographically compare w/ tolerance
        {
            value_type const & x = _lhs.p.x + eps_;
            value_type const & y = _lhs.p.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.p.x, _rhs.p.y);
        }

    };

    using vertices = std::set< vertex, vertex_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::list< edge >;
    using pedge = typename edges::iterator;

    using cell = std::deque< pedge >; // front() - left, back() - right
    using cells = std::map< site, cell >;
    using pcell = typename cells::iterator;

    vertex_less const vertex_less_{eps};
    vertices vertices_{vertex_less_};
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_;

private :

    struct endpoint
    {

        pcell l, r;
        pedge e;

    };

    value_type directrix = std::numeric_limits< value_type >::quiet_NaN();

    struct local_endpoints
    {

        endpoint const * l = nullptr;
        endpoint const * m = nullptr;
        endpoint const * r = nullptr;

        void clear() { l = m = r = nullptr; }

    };

    struct endpoint_less
    {

        value_type const & eps_;

        value_type const & directrix_;

        value_type // ordinate
        intersect(point_type const & u,
                  point_type const & v,
                  value_type _directrix) const
        {
            {
                bool const degenerated_ = !(v.x + eps_ < _directrix);
                if (!(u.x + eps_ < _directrix)) {
                    if (degenerated_) {
                        assert(u.y + eps_ < v.y); // u != v
                        return (u.y + v.y) / value_type(2);
                    } else {
                        return u.y;
                    }
                } else if (degenerated_) {
                    return v.y;
                }
            }
            value_type ud = u.x - _directrix;
            value_type vd = v.x - _directrix;
            value_type ub = u.y / ud; // -b
            value_type vb = v.y / vd; // -b
            ud += ud;
            vd += vd;
            _directrix *= _directrix;
            auto const calc_c = [&] (auto const & w, value_type const & wd)
            {
                return (w.x * w.x + w.y * w.y - _directrix) / wd;
            };
            value_type uc = calc_c(u, ud);
            value_type vc = calc_c(v, vd);
            value_type b = vb - ub; // -b
            value_type c = vc - uc;
            if ((u.x + eps_ < v.x) || (v.x + eps_ < u.x)) {
                value_type a = (ud - vd) / (ud * vd);
                a += a;
                value_type D = b * b - (a + a) * c;
                assert(!(D < value_type(0)));
                using std::sqrt;
                return (b + sqrt(std::move(D))) / a;
            } else { // a ~= 0
                return c / b; // -c / b
            }
        }

        local_endpoints const & lep_;

        bool operator () (endpoint const & _lhs, endpoint const & _rhs) const
        {
            if (&_lhs == &_rhs) {
                return false;
            }
            if (lep_.m) {
                if (lep_.m == &_lhs) {
                    if (lep_.l == &_rhs) {
                        return false;
                    }
                    if (!lep_.r) {
                        return false;
                    }
                    if (lep_.r == &_rhs) {
                        return true;
                    }
                } else if (lep_.m == &_rhs) {
                    if (lep_.r == &_lhs) {
                        return false;
                    }
                    if (!lep_.l) {
                        return false;
                    }
                    if (lep_.l == &_rhs) {
                        return true;
                    }
                }
            }
            if (_lhs.r == _rhs.l) {
                return true;
            }
            if (_lhs.l == _rhs.r) {
                return false;
            }
            return intersect(*_lhs.l->first, *_lhs.r->first, directrix_) + eps_ < intersect(*_rhs.l->first, *_rhs.r->first, directrix_);
        }

    };

    using endpoints = std::map< endpoint, pvertex, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    local_endpoints lep;
    endpoint_less const endpoint_less_{eps, directrix, lep};
    endpoints endpoints_{endpoint_less_};

    struct event_less
    {

        value_type const & eps_;

        bool operator () (pvertex const & _lhs, pvertex const & _rhs) const
        {
            vertex const & lhs_ = *_lhs;
            value_type const & x = lhs_.x() + eps_;
            value_type const & y = lhs_.y() + eps_;
            vertex const & rhs_ = *_rhs;
            value_type const & xx = rhs_.x();
            return std::tie(x, y) < std::tie(xx, rhs_.y());
        }

    };

    using events = std::multimap< pvertex, pendpoint const, event_less >;
    using pevent = typename events::const_iterator;

    event_less const event_less_{eps};
    events events_{event_less_};

    void
    begin_cell(site const &, cell &)
    {

    }

    std::pair< pvertex, bool >
    make_vertex(point_type const & u,
                point_type const & v,
                point_type const & w)
    {
        value_type A = v.x - u.x;
        value_type B = v.y - u.y;
        value_type C = w.x - u.x;
        value_type D = w.y - u.y;
        value_type G = B * (w.x - v.x) - A * (w.y - v.y);
        if (!(eps * eps < G)) {
            // 1.) non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return {nov, false};
        }
        G += G;
        value_type E = A * (u.x + v.x) + B * (u.y + v.y);
        value_type F = C * (u.x + w.x) + D * (u.y + w.y);
        // x, y - circumcenter
        using std::sqrt;
        auto const norm = [&] (auto const & ll, auto const & rr) -> value_type
        {
            value_type dx = rr.x - ll.x;
            value_type dy = rr.y - ll.y;
            return sqrt(dx * dx + dy * dy);
        };
        value_type e = norm(u, v);
        value_type f = norm(v, w);
        value_type g = norm(w, u);
        value_type V = (e + f - g) * (e + g - f) * (f + g - e);
        assert(eps < V); // triangle inequality
        return vertices_.insert({{(B * F - D * E) / G, (C * E - A * F) / G},
                                 (e * f * g) / sqrt(std::move(V) * (e + f + g))});
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & a = *l;
        auto & b = *r;
        assert(a.first.r == b.first.l);
        pvertex const v = make_vertex(*a.first.l->first,
                                      *a.first.r->first,
                                      *b.first.r->first).first;
        if (v != nov) {
            a.second = v; // disable corresponding event if a.second != nov
            b.second = v; // disable corresponding event
            events_.insert({v, l});
            events_.insert({v, r});
        }
    }

    void
    finish_edge(edge & _edge, pvertex const v) const
    {
        assert(v != nov);
        assert(_edge.e != v);
        if (_edge.b == nov) {
            _edge.b = v;
            auto const & l = *_edge.l;
            auto const & r = *_edge.r;
            point_type const & p = v->p;
            if (r.x < l.x) {
                if (p.y < l.y) {
                    return;
                }
            } else if (l.x < r.x) {
                if (r.y < p.y) {
                    return;
                }
            } else {
                if (r.y < l.y) {
                    return; // unrecheable
                }
            }
            std::swap(_edge.l, _edge.r);
        } else {
            assert(_edge.b != v);
            assert(_edge.e == nov);
            _edge.e = v;
        }
    }

    void
    collapse(pendpoint l, pendpoint r, pvertex const v)
    {
        pcell const lcell = l->first.l;
        pcell const rcell = r->first.r;
        ++r;
        for (pendpoint ep = l; ep != r; ++ep) {
            finish_edge(*ep->first.e, v);
        }
        endpoints_.erase(l, r);
        l = endpoints_.insert(r, {{lcell, rcell, edges_.insert(std::cend(edges_), {lcell->first, rcell->first, v, nov})}, nov});
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != std::end(endpoints_)) {
            check_event(l, r);
        }
    }

    void
    finish_edges(pevent const l, pevent const r, pvertex const v)
    {
        assert(l != r);
        assert(std::next(l) != r);
        auto const pendpoint_less_ = [&] (auto const & _lhs, auto const & _rhs)
        {
            assert(_lhs.first == v);
            assert(_rhs.first == v);
            assert(_lhs.second->second == v);
            assert(_rhs.second->second == v);
            return endpoint_less_(_lhs.second->first, _rhs.second->first);
        };
        auto const lr = std::minmax_element(l, r, pendpoint_less_);
        assert(lr.first != lr.second);
        events_.erase(l, r);
        collapse(lr.first->second, lr.second->second, v);
    }

    bool
    prior(vertex const & _lhs, point_type const & _rhs) const
    {
        value_type const & x = _lhs.x() + eps;
        value_type const & y = _lhs.y() + eps;
        return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
    }

public :

    bool
    operator () (point_iterator beg, point_iterator const end)
    {
        if (beg == end) {
            return true;
        }
        do {
            if (!cells_.insert({beg, {}}).second) {
                return false;
            }
        } while (++beg != end);
        for (auto & cell_ : cells_) {
            while (!events_.empty()) {
                pevent const e = std::cbegin(events_);
                pvertex const v = e->first;
                if (!prior(*v, *cell_.first)) {
                    break;
                }
                finish_edges(e, events_.upper_bound(v), v);
            }
            begin_cell(cell_.first, cell_.second);
        }
        while (!events_.empty()) {
            pevent const e = std::cbegin(events_);
            pvertex const v = e->first;
            finish_edges(e, events_.upper_bound(v), v);
        }
        return true;
    }

};
