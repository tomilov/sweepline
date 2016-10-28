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

    vertices vertices_{point_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;

private :

    struct endpoint
    {

        point_iterator l, r;
        pedge e;

        bool operator == (endpoint const & ep) const
        {
            if (this == &ep) {
                return true;
            }
            return (l == ep.l) && (r == ep.r) && (e == ep.e);
        }

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
        intersect(endpoint const & ep, value_type _directrix) const
        {
            return intersect(*ep.l, *ep.r, std::move(_directrix));
        }

        bool operator () (endpoint const & l, endpoint const & r) const
        {
            // during sweepline motion arcs shrinks and growz, but relative y-position of endpoints remains the same
            // endpoints removed strictly before violation of this invariant to prevent its occurrence
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            if (l == r) {
                return false;
            }
            if (l.l == r.l) {
                return true;
            }
            if (r.r == l.r) {
                return true;
            }/*
            if (std::max(*l.l, *l.r).y < std::max(*r.l, *r.r).y) {
                return true;
            } else {
                return false;
            }*/
            throw std::logic_error{""};
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
        //value_type R = circumradius(hypot(A, B), hypot(C, D), hypot(E, F));
        value_type R = hypot(x - a.x, y - a.y);
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
        pvertex const v = make_vertex(*ll.first.l, *ll.first.r, *rr.first.r);
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
    trunc_edge(edge & ee, pvertex const v) const
    {
        assert(v != nov);
        if (ee.b == nov) {
            if (ee.e == nov) { // orientate if needed:
                point_type const & l = *ee.l;
                point_type const & r = *ee.r;
                point_type const & p = v->p;
                if (r.x < l.x) {
                    if (p.y < l.y) {
                        ee.b = v;
                        return;
                    }
                } else if (l.x < r.x) {
                    if (r.y < p.y) {
                        ee.b = v;
                        return;
                    }
                } else {
                    assert(!(r.y < l.y));
                }
                ee.e = v;
            } else {
                assert(ee.e != v);
                ee.b = v;
            }
        } else {
            assert(ee.b != v);
            assert(ee.e == nov);
            ee.e = v;
        }
    }

    pedge
    add_edge(point_iterator const l, point_iterator const r)
    {
        return edges_.insert(std::cend(edges_), {l, r, nov, nov});
    }

    void
    finish_endpoints(pendpoint l, pendpoint const r, pvertex const v)
    {
        point_iterator const lc = l->first.l;
        point_iterator const rc = std::prev(r)->first.r;
        for (pendpoint ep = l; ep != r; ++ep) {
            assert(ep->second == v);
            trunc_edge(*ep->first.e, v);
        }
        endpoints_.erase(l, r);
        pedge const e = add_edge(lc, rc);
        trunc_edge(*e, v);
        l = endpoints_.insert(r, {{lc, rc, e}, nov});
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != noep) {
            check_event(l, r);
        }
    }

    void
    make_first_edge(point_iterator const l, point_iterator const r)
    {
        assert(endpoints_.empty());
        pedge const ee = add_edge(l, r);
        pendpoint const le = endpoints_.insert(noep, {{l, r, ee}, nov});
        pendpoint const re = endpoints_.insert(noep, {{r, l, ee}, nov});
        assert(std::next(le) == re);
    }

    void
    begin_cell(point_iterator const p)
    {
        auto lr = endpoints_.equal_range(*p);
        if (lr.first == lr.second) {
            if (lr.first == noep) {
                lr.first = std::prev(noep);
                point_iterator const l = lr.first->first.r;
                pedge const e = add_edge(l, p);
                lr.second = endpoints_.insert(noep, {{p, l, e}, nov});
                endpoints_.insert(noep, {{l, p, e}, nov});
                check_event(lr.first, lr.second);
            } else if (lr.first == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                point_iterator const r = lr.second->first.l;
                pedge const e = add_edge(p, r);
                endpoints_.insert(lr.second, {{r, p, e}, nov});
                lr.first = endpoints_.insert(lr.second, {{p, r, e}, nov});
                check_event(lr.first, lr.second);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --lr.first;
                point_iterator const c = lr.first->first.r;
                assert(c == lr.second->first.l);
                pedge const e = add_edge(c, p);
                pendpoint const l = endpoints_.insert(lr.second, {{c, p, e}, nov});
                pendpoint const r = endpoints_.insert(lr.second, {{p, c, e}, nov});
                assert(std::next(l) == r);
                check_event(lr.first, l);
                check_event(r, lr.second);
            }
        } else { // many arc collapsing right here, event (equivalent to the current site p) coming on the next step
            pvertex const v = lr.first->second;
            events_.erase(v);
            finish_endpoints(lr.first, lr.second, v);
        }
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
    operator () (point_iterator l, point_iterator const r)
    {
        assert(std::is_sorted(l, r, point_less{eps}));
        if (l == r) {
            return;
        }
        if (std::next(l) == r) {
            return;
        }
        make_first_edge(l, ++l);
        while (++l != r) {
            while (!events_.empty()) {
                pevent const e = std::cbegin(events_);
                if (!prior(*e->first, *l)) {
                    break;
                }
                finish_cells(e);
            }
            begin_cell(l);
        }
        while (!events_.empty()) {
            finish_cells(std::cbegin(events_));
        }
        endpoints_.clear();
    }

};
