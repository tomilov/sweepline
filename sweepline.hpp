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

    static_assert(std::is_same< decltype(std::declval< point_type >().x), decltype(std::declval< point_type >().y) >::value,
                  "point_type format error");

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct site_less
    {

        bool operator () (point_iterator const l, point_iterator const r) const
        {
            return *l < *r;
        }

    };

    struct vertex // denote circumscribed circle
    {

        point_type p; // circumcenter
        value_type R; // radius

        value_type x() const { return p.x + R; }
        value_type const & y() const { return p.y; }

    };

    struct vertex_less // lexicographically compare w/ tolerance
    {

        value_type const & eps_;

        bool operator () (vertex const & l, vertex const & r) const
        {
            return operator () (l.p, r.p);
        }

        bool operator () (point_type const & l, point_type const & r) const
        {
            value_type const & x = l.x + eps_;
            value_type const & y = l.y + eps_;
            return std::tie(x, y) < std::tie(r.x, r.y);
        }

    };

    using vertices = std::set< vertex, vertex_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        point_iterator l, r;
        pvertex b, e;

    };

    using edges = std::list< edge >;
    using pedge = typename edges::iterator;

    using cells = std::map< point_iterator, std::deque< pedge >, site_less >;
    using pcell = typename cells::iterator;

    vertices vertices_{vertex_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_;

private :

    struct endpoint
    {

        pcell l, r;
        pedge e;

    };

    struct local_insert_hint
    {

        endpoint const * l = nullptr;
        endpoint const * m = nullptr;
        endpoint const * r = nullptr;

        void clear() { l = m = r = nullptr; }

    };

    struct endpoint_less
    {

        value_type const & eps_;

        value_type
        intersect(point_type const & l,
                  point_type const & r,
                  value_type && _directrix) const
        {
            {
                bool const degenerated_ = !(r.x + eps_ < _directrix);
                if (!(l.x + eps_ < _directrix)) {
                    if (degenerated_) {
                        assert(l.y + eps_ < r.y); // l != r
                        return (l.y + r.y) / value_type(2);
                    } else {
                        return l.y;
                    }
                } else if (degenerated_) {
                    return r.y;
                }
            }
            value_type ld = l.x - _directrix;
            value_type rd = r.x - _directrix;
            value_type lb = l.y / ld; // -b
            value_type rb = r.y / rd; // -b
            ld += ld;
            rd += rd;
            _directrix *= _directrix;
            auto const calc_c = [&] (point_type const & p, value_type const & d)
            {
                return (p.x * p.x + p.y * p.y - _directrix) / d;
            };
            value_type lc = calc_c(l, ld);
            value_type rc = calc_c(r, rd);
            value_type b = rb - lb; // -b
            value_type c = rc - lc;
            if ((l.x + eps_ < r.x) || (r.x + eps_ < l.x)) {
                value_type a = (ld - rd) / (ld * rd);
                a += a;
                value_type D = b * b - (a + a) * c;
                assert(!(D < value_type(0)));
                using std::sqrt;
                return (b + sqrt(D)) / a;
            } else { // a ~= 0
                return c / b; // -c / b
            }
        }

        value_type
        intersect(endpoint const & lr, value_type _directrix) const
        {
            return intersect(*lr.l->first, *lr.r->first, std::move(_directrix));
        }

        local_insert_hint const & hint_;

        bool operator () (endpoint const & l, endpoint const & r) const
        {
            if (&l == &r) {
                return false;
            }
            if (hint_.m) {
                if (hint_.m == &l) {
                    if (!hint_.r) {
                        return false;
                    }
                    if (hint_.l == &r) {
                        return false;
                    }
                    if (hint_.r == &r) {
                        return true;
                    }
                } else if (hint_.m == &r) {
                    if (!hint_.l) {
                        return false;
                    }
                    if (hint_.r == &l) {
                        return false;
                    }
                    if (hint_.l == &r) {
                        return true;
                    }
                }
            }
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            return std::max(*l.l->first, *l.r->first, vertex_less{eps_}).y < std::max(*r.l->first, *r.r->first, vertex_less{eps_}).y;
        }

        using is_transparent = void;

        bool operator () (vertex const & l, endpoint const & r) const
        {
            return l.y() + eps_ < intersect(r, l.x());
        }

        bool operator () (endpoint const & l, vertex const & r) const
        {
            return intersect(l, r.x()) + eps_ < r.y();
        }

        bool operator () (point_type const & l, endpoint const & r) const
        {
            return l.y + eps_ < intersect(r, l.x);
        }

        bool operator () (endpoint const & l, point_type const & r) const
        {
            return intersect(l, r.x) + eps_ < r.y;
        }

    };

    using endpoints = std::map< endpoint, pvertex, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    local_insert_hint hint;
    endpoint_less const endpoint_less_{eps, hint};
    endpoints endpoints_{endpoint_less_};
    pendpoint const noep = std::end(endpoints_);

    struct event_less
    {

        value_type const & eps_;

        bool operator () (pvertex const & l, pvertex const & r) const
        { // uncomment if not all the events in equal range are for equally the same vertex
            vertex const & ll = *l;
            value_type const & x = ll.x() + eps_;
            value_type const & y = ll.y() + eps_;
            //value_type const & px = lhs_.p.x + eps_;
            vertex const & rr = *r;
            value_type const & vx = rr.x();
            return std::tie(x, y/*, px*/) < std::tie(vx, rr.y()/*, rhs_.p.x*/);
        }

    };

    using events = std::map< pvertex, pendpoint, event_less >;
    using pevent = typename events::const_iterator;

    events events_{event_less{eps}};

    void
    begin_cell(point_iterator const p)
    {
        assert(cells_.find(p) == cells_.end());
        pcell const c = cells_.insert({p, {}}).first;
        (void)c; // TODO: implement
    }

    value_type
    circumradius(value_type a,
                 value_type b,
                 value_type c) const
    {
        value_type V = (a + b - c) * (a + c - b) * (b + c - a);
        assert(eps < V); // triangle inequality
        using std::sqrt;
        return (a * b * c) / sqrt(V * (a + b + c));
    }

    pvertex
    make_vertex(point_type const & u,
                point_type const & v,
                point_type const & w)
    {
        value_type A = v.x - u.x;
        value_type B = v.y - u.y;
        value_type C = w.x - v.x;
        value_type D = w.y - v.y;
        value_type G = B * C - A * D;
        if (!(eps * eps < G)) {
            // 1.) non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return nov;
        }
        value_type E = w.x - u.x;
        value_type F = w.y - u.y;
        value_type M = A * (u.x + v.x) + B * (u.y + v.y);
        value_type N = E * (u.x + w.x) + F * (u.y + w.y);
        G += G;
        // circumcenter:
        value_type x = (B * N - F * M) / G;
        value_type y = (E * M - A * N) / G;
        using std::hypot;
        value_type R = circumradius(hypot(A, B), hypot(C, D), hypot(E, F));
        return vertices_.insert({{x, y}, R}).first;
    }

    std::pair< pendpoint, pendpoint >
    endpoint_range(pendpoint const ep, pvertex const v)
    {
        if (ep == noep) {
            return endpoints_.equal_range(*v);
        } else {
            assert(ep != std::begin(endpoints_));
            assert(ep != noep);
            assert(std::prev(ep)->second == v);
            assert(ep->second == v);
            return {std::prev(ep), std::next(ep)};
        }
    }

    void
    delete_event(pvertex const v)
    {
        assert(events_.find(v) != std::end(events_));
        auto const e = events_.find(v);
        auto const lr = endpoint_range(e->second, v);
        events_.erase(e);
        for (auto ep = lr.first; ep != lr.second; ++ep) {
            ep->second = nov;
        }
        vertices_.erase(v);
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l); // small redundancy (x2 memory to point cells from beach line) - necessary evil
        pvertex const v = make_vertex(*ll.first.l->first, *ll.first.r->first, *rr.first.r->first);
        if (v != nov) {
            if (ll.second != nov) {
                assert(rr.second == nov);
                if (event_less{eps}(v, ll.second)) {
                    delete_event(ll.second);
                } else {
                    vertices_.erase(v);
                    return;
                }
            } else if (rr.second != nov) {
                assert(ll.second == nov);
                if (event_less{eps}(v, rr.second)) {
                    delete_event(rr.second);
                } else {
                    vertices_.erase(v);
                    return;
                }
            }
            assert(ll.second == nov);
            assert(rr.second == nov);
            ll.second = rr.second = v;
            auto const e = events_.insert({v, r});
            if (!e.second) {
                e.first->second = noep;
            }
        }
    }

    void
    finish_edge(edge & _edge, pvertex const v) const
    {
        assert(v != nov);
        assert(_edge.e == nov);
        if (_edge.b == nov) {
            _edge.b = v;
            // adjust orinetation if needed:
            point_type const & l = *_edge.l;
            point_type const & r = *_edge.r;
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
                assert(!(r.y < l.y));
            }
            std::swap(_edge.l, _edge.r);
        } else {
            assert(_edge.b != v);
            _edge.e = v;
        }
    }

    void
    finish_edges(pevent const e)
    {
        pvertex const v = e->first;
        auto lr = endpoint_range(e->second, v);
        events_.erase(e);
        pcell const lcell = lr.first->first.l;
        pcell const rcell = std::prev(lr.second)->first.r;
        for (pendpoint ep = lr.first; ep != lr.second; ++ep) {
            finish_edge(*ep->first.e, v);
        }
        endpoints_.erase(lr.first, lr.second);
        pedge const edge_ = edges_.insert(std::cend(edges_), {lcell->first, rcell->first, v, nov});
        lr.first = endpoints_.insert(lr.second, {{lcell, rcell, edge_}, nov});
        lcell->second.push_front(edge_); // ccw
        rcell->second.push_back(edge_);
        if (lr.first != std::begin(endpoints_)) {
            check_event(std::prev(lr.first), lr.first);
        }
        if (lr.second != noep) {
            check_event(lr.first, lr.second);
        }
    }

    bool
    prior(vertex const & l, point_type const & r) const
    {
        value_type const & x = l.x() + eps;
        value_type const & y = l.y() + eps;
        return std::tie(x, y) < std::tie(r.x, r.y);
    }

public :

    void
    operator () (point_iterator const beg, point_iterator const end)
    {
        for (auto p = beg; p != end; ++p) {
            while (!events_.empty()) {
                pevent const e = std::cbegin(events_);
                if (!prior(*e->first, *p)) {
                    break;
                }
                finish_edges(e);
            }
            begin_cell(p);
        }
        while (!events_.empty()) {
            finish_edges(std::cbegin(events_));
        }
    }

};
