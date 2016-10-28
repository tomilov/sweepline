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

    struct vertex // denote circumscribed circle
    {

        point_type p; // circumcenter
        value_type R; // radius

        value_type x() const { return p.x + R; }
        value_type const & y() const { return p.y; }

    };

    struct point_less // lexicographically compare w/ tolerance
    {

        value_type const & eps_;

        bool operator () (point_type const & l, point_type const & r) const
        {
            value_type const & x = l.x + eps_;
            value_type const & y = l.y + eps_;
            return std::tie(x, y) < std::tie(r.x, r.y);
        }

        bool operator () (point_iterator const l, point_iterator const r) const
        {
            return operator () (*l, *r);
        }

        bool operator () (vertex const & l, vertex const & r) const
        {
            return operator () (l.p, r.p);
        }

    };

    using vertices = std::set< vertex, point_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        point_iterator l, r;
        pvertex b, e;

    };

    using edges = std::list< edge >;
    using pedge = typename edges::iterator;

    using cells = std::map< point_iterator, std::deque< pedge >, point_less >;
    using pcell = typename cells::iterator;

    vertices vertices_{point_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_{point_less{eps}};

private :

    struct endpoint
    {

        pcell l, r;
        pedge e;

        bool operator == (endpoint const & ep) const
        {
            if (this == &ep) {
                return true;
            }
            return (l == ep.l) && (r == ep.r) && (assert(e == ep.e), true);
        }

    };

    struct endpoint_less
    {

        value_type const & eps_;

        bool operator () (endpoint const & l, endpoint const & r) const
        {
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            if (l == r) {
                return false;
            }/*
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
            }*/
            // during sweepline motion arcs shrinks and growz, but relative y-position of endpoints remains the same
            // endpoints removed strictly before violation of this invariant to prevent its occurrence
            return std::max(*l.l->first, *l.r->first, point_less{eps_}).y < std::max(*r.l->first, *r.r->first, point_less{eps_}).y;
        }

        using is_transparent = void;

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

    endpoints endpoints_{endpoint_less{eps}};
    pendpoint const noep = std::end(endpoints_);

    struct event_less
    {

        value_type const & eps_;

        bool operator () (pvertex const & l, pvertex const & r) const
        {
            vertex const & ll = *l;
            value_type const & lx = ll.x() + eps_;
            value_type const & ly = ll.y() + eps_;
            vertex const & rr = *r;
            value_type const & rx = rr.x();
            return std::tie(lx, ly) < std::tie(rx, rr.y());
        }

    };

    using events = std::map< pvertex, pendpoint, event_less >;
    using pevent = typename events::const_iterator;

    events events_{event_less{eps}};

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
    make_vertex(point_type const & a,
                point_type const & b,
                point_type const & c)
    {
        value_type A = b.x - a.x;
        value_type B = b.y - a.y;
        value_type C = c.x - b.x;
        value_type D = c.y - b.y;
        value_type G = B * C - A * D;
        if (!(eps * eps < G)) {
            // 1.) G is negative: non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return nov;
        }
        value_type E = c.x - a.x;
        value_type F = c.y - a.y;
        value_type M = A * (a.x + b.x) + B * (a.y + b.y);
        value_type N = E * (a.x + c.x) + F * (a.y + c.y);
        G += G;
        // circumcenter:
        value_type x = (B * N - F * M) / G;
        value_type y = (E * M - A * N) / G;
        using std::hypot;
        value_type R = circumradius(hypot(A, B), hypot(C, D), hypot(E, F));
        //value_type R = hypot(x - a.x, y - a.y);
        return vertices_.insert({{x, y}, R}).first;
    }

    std::pair< pendpoint, pendpoint >
    endpoint_range(pendpoint const ep, pvertex const v)
    {
        if (ep == noep) {
            return endpoints_.equal_range(*v);
        }
        assert(ep != std::begin(endpoints_));
        assert(ep != noep);
        assert(std::prev(ep)->second == v);
        assert(ep->second == v);
        return {std::prev(ep), std::next(ep)};
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

    pvertex
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
                    return nov;
                }
            } else if (rr.second != nov) {
                assert(ll.second == nov);
                if (event_less{eps}(v, rr.second)) {
                    delete_event(rr.second);
                } else {
                    vertices_.erase(v);
                    return nov;
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
        return v;
    }

    void
    trunc_edge(edge & _edge, pvertex const v) const
    {
        assert(v != nov);
        if (_edge.b == nov) {
            if (_edge.e == nov) { // orientate if needed:
                point_type const & l = *_edge.l;
                point_type const & r = *_edge.r;
                point_type const & p = v->p;
                if (r.x < l.x) {
                    if (p.y < l.y) {
                        _edge.b = v;
                        return;
                    }
                } else if (l.x < r.x) {
                    if (r.y < p.y) {
                        _edge.b = v;
                        return;
                    }
                } else {
                    assert(!(r.y < l.y));
                }
                _edge.e = v;
            } else {
                assert(_edge.e != v);
                _edge.b = v;
            }
        } else {
            assert(_edge.b != v);
            assert(_edge.e == nov);
            _edge.e = v;
        }
    }

    void
    finish_endpoints(pendpoint l, pendpoint const r, pvertex const v)
    {
        pcell const lc = l->first.l;
        pcell const rc = std::prev(r)->first.r;
        for (pendpoint ep = l; ep != r; ++ep) {
            assert(ep->second == v);
            trunc_edge(*ep->first.e, v);
        }
        endpoints_.erase(l, r);
        pedge const edge_ = edges_.insert(std::cend(edges_), {lc->first, rc->first, v, nov});
        l = endpoints_.insert(r, {{lc, rc, edge_}, nov});
        lc->second.push_front(edge_); // ccw
        rc->second.push_back(edge_);
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != noep) {
            check_event(l, r);
        }
    }

    pcell
    begin_cell(point_iterator const p)
    {
        pcell const c = cells_.insert(std::cend(cells_), {p, {}});
        auto lr = endpoints_.equal_range(*p);
        if (lr.first == noep) { // begginning
            if (endpoints_.empty()) {
                pcell const lc = std::begin(cells_);
                if (c == lc) { // first site
                    assert(cells_.size() == 1);
                } else { // second site
                    assert(cells_.size() == 2);
                    assert(std::next(lc) == c);
                    auto & lcell = *lc;
                    pedge const e = edges_.insert(std::cend(edges_), {lcell.first, p, nov, nov});
                    endpoints_.insert({{lc, c, e}, nov});
                    lcell.second.push_front(e);
                    c->second.push_back(e);
                }
            } else { // appending to the rightmost endpoint
                assert(2 < cells_.size()); // another site
                lr.first = std::prev(noep);
                pcell const lc = lr.first->first.r;
                auto & lcell = *lc;
                pedge const le = edges_.insert(std::cend(edges_), {lcell.first, p, nov, nov});
                lr.second = endpoints_.insert(noep, {{lc, c, le}, nov});
                lcell.second.push_front(le);
                c->second.push_back(le);
                if (lr.first != std::begin(endpoints_)) {
                    check_event(std::prev(lr.first), lr.first);
                }
                check_event(lr.first, lr.second);
            }
        } else {
            if (lr.first == lr.second) { // add a single endpoint
                if (lr.first == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                    pcell const rc = lr.second->first.l;
                    auto & rcell = *rc;
                    pedge const re = edges_.insert(std::cend(edges_), {p, rcell.first, nov, nov});
                    lr.first = endpoints_.insert(lr.second, {{c, rc, re}, nov});
                    c->second.push_front(re);
                    rcell.second.push_back(re);
                    check_event(lr.first, lr.second);
                    if (2 < endpoints_.size()) {
                        check_event(lr.second, std::next(lr.second));
                    }
                } else { // insert in the middle of the beachline (hottest branch in general case)
                    lr.first = std::prev(lr.second);
                    pcell const lc = lr.first->first.r;
                    pcell const rc = lr.second->first.l;
                    auto & lcell = *lc;
                    auto & mcell = *c;
                    auto & rcell = *rc;
                    pedge const le = edges_.insert(std::cend(edges_), {lcell.first, p, nov, nov});
                    pedge const re = edges_.insert(std::cend(edges_), {p, rcell.first, nov, nov});
                    pendpoint const l = endpoints_.insert(lr.second, {{lc, c, le}, nov});
                    pendpoint const r = endpoints_.insert(lr.second, {{c, rc, re}, nov});
                    lcell.second.push_front(le);
                    mcell.second.push_back(le);
                    mcell.second.push_front(re);
                    rcell.second.push_back(re);
                    check_event(lr.first, l);
                    check_event(r, lr.second);
                }
            } else { // many arc collapsing right here, event (equivalent to the current site) coming on the next step
                pvertex const v = lr.first->second;
                events_.erase(v);
                finish_endpoints(lr.first, lr.second, v);
            }
        }
        return c;
    }

    void
    finish_cells(pevent const e)
    {
        pvertex const v = e->first;
        auto lr = endpoint_range(e->second, v);
        events_.erase(e);
        finish_endpoints(lr.first, lr.second, v);
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
    operator () (point_iterator const l, point_iterator const r)
    {
        assert(std::is_sorted(l, r, point_less{eps}));
        for (auto p = l; p != r; ++p) {
            while (!events_.empty()) {
                pevent const e = std::cbegin(events_);
                if (!prior(*e->first, *p)) {
                    break;
                }
                finish_cells(e);
            }
            begin_cell(p);
        }
        while (!events_.empty()) {
            finish_cells(std::cbegin(events_));
        }
    }

};
