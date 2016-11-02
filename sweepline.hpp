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
#include <iomanip>
#include <ostream>
#include <iostream>

#include <cassert>
#include <cmath>

template< typename site,
          typename point_less,
          typename point,
          typename value_type >
struct sweepline
{

    static_assert(std::is_same< decltype(std::declval< point >().x), decltype(std::declval< point >().y) >::value,
                  "point format error");

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct vertex // denote circumscribed circle
    {

        point c; // circumcenter
        value_type R; // circumradius

        value_type x() const { return c.x + R; }
        value_type const & y() const { return c.y; }

        friend
        std::ostream &
        operator << (std::ostream & _out, vertex const & v)
        {
            return _out << "v{" << v.c << ", " << v.R << '}';
        }

    };

    struct vertex_less
    {

        value_type const & eps_;

        bool operator () (vertex const & l, vertex const & r) const
        {
            return point_less{eps_}(l.c, r.c);
        }

    };

    using vertices = std::set< vertex, vertex_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        site l, r;
        pvertex b, e;

        friend
        std::ostream &
        operator << (std::ostream & _out, edge const & ed)
        {
            return _out << "e{" << *ed.l << ", " << *ed.r << ", " << *ed.b << ", " << *ed.e << '}';
        }

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::iterator;

    vertices vertices_{vertex_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;

private :

    struct endpoint
    {

        site l, r;
        pedge e;

        bool operator == (endpoint const & ep) const
        {
            if (this == &ep) {
                return true;
            }
            if ((l == ep.l) && (r == ep.r)) {
                assert(e == ep.e);
                return true;
            } else {
                return false;
            }
        }

        friend
        std::ostream &
        operator << (std::ostream & _out, endpoint const & ep)
        {
            return _out << "ep{" << *ep.l << ", " << *ep.r << ", " << *ep.e << '}' << std::endl;
        }

    };

    static
    value_type
    intersect(point const & p,
              value_type const & y,
              value_type const & directrix)
    {
        assert(p.x < directrix);
        value_type d = p.x - directrix;
        return (y * (y - (p.y + p.y)) + (p.x * p.x + p.y * p.y - directrix * directrix)) / (d + d);
    }

    struct endpoint_less
    {

        value_type const & eps_;

        value_type
        intersect(point const & l,
                  point const & r,
                  value_type const & directrix) const
        {
            {
                bool const rdegenerated = !(r.x + eps_ < directrix);
                if (!(l.x + eps_ < directrix)) {
                    if (rdegenerated) {
                        assert((l.y + eps_ < r.y) || (r.y + eps_ < l.y)); // l != r
                        return (l.y + r.y) / value_type(2);
                    } else {
                        return l.y;
                    }
                } else if (rdegenerated) {
                    return r.y;
                }
            }
            value_type ld = l.x - directrix;
            value_type rd = r.x - directrix;
            value_type lb = l.y / ld; // -b
            value_type rb = r.y / rd; // -b
            ld += ld;
            rd += rd;
            auto const calc_c = [&] (point const & p, value_type const & d)
            {
                return (p.x * p.x + p.y * p.y - directrix * directrix) / d;
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
                value_type y = (b + sqrt(D)) / a;
                {
                    /*
                    value_type x0 = (l.x + r.x) / value_type(2);
                    value_type y0 = (l.y + r.y) / value_type(2);
                    value_type dx = l.y - r.y;
                    value_type dy = r.x - l.x;
                    value_type xx = x0 + (y - y0) * dx / dy;
                    */
                    value_type xa = sweepline::intersect(l, y, directrix);
                    value_type xb = sweepline::intersect(r, y, directrix);
                    value_type da = std::hypot(xa - l.x, y - l.y);
                    value_type db = std::hypot(xb - r.x, y - r.y);
                    asm volatile ("nop");
                }
                return y;
            } else { // a ~= 0
                return c / b; // -c / b
            }
        }

        value_type
        intersect(endpoint const & ep, value_type const & directrix) const
        {
            return intersect(*ep.l, *ep.r, directrix);
        }

        bool operator () (endpoint const & l, endpoint const & r) const
        {
            // during sweepline motion arcs shrinks and growz, but relative y-position of endpoints remains the same
            // endpoints removed strictly before violation of this invariant to prevent its occurrence
            assert(!(l == r));
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            if (l.l == r.l) {
                return true;
            }
            if (r.r == l.r) {
                return true;
            }
            throw std::logic_error{"undefined behaviour"};
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

        void print_delta(vertex const & v, endpoint const & ep) const
        {
            std::cerr << v.y() - intersect(ep, v.x()) << std::endl;
        }

        bool operator () (point const & l, endpoint const & r) const
        {
            return l.y + eps_ < intersect(r, l.x);
        }

        bool operator () (endpoint const & l, point const & r) const
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

        bool operator () (vertex const & l, vertex const & r) const
        {
            value_type const & lx = l.x() + eps_;
            value_type const & ly = l.y() + eps_;
            value_type const & rx = r.x();
            return std::tie(lx, ly) < std::tie(rx, r.y());
        }

        bool operator () (pvertex const l, pvertex const r) const
        {
            if (l == r) {
                return false;
            }
            return operator () (*l, *r);
        }

    };

    using events = std::multimap< pvertex, pendpoint, event_less >;
    using pevent = typename events::iterator;

    events events_{event_less{eps}};
    pevent const noe = std::end(events_);

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

    std::pair< pvertex, bool >
    make_vertex(point const & a,
                point const & b,
                point const & c)
    {
        value_type A = b.x - a.x;
        value_type B = b.y - a.y;
        value_type C = c.x - b.x;
        value_type D = c.y - b.y;
        value_type G = B * C - A * D;
        if (!(eps * eps < G)) {
            // 1.) G is negative: non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return {nov, false};
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
#if 0
        value_type R = circumradius(hypot(A, B), hypot(C, D), hypot(E, F));
#elif 0
        value_type R = (hypot(x - a.x, y - a.y) + hypot(x - b.x, y - b.y) + hypot(x - c.x, y - c.y)) / value_type(3);
#else
        value_type R = hypot(x - a.x, y - a.y);
#endif
        {
            value_type x1 = intersect(a, y, x + R);
            value_type x2 = intersect(b, y, x + R);
            value_type x3 = intersect(c, y, x + R);
            value_type y1 = endpoint_less{eps}.intersect(a, b, x + R);
            value_type y2 = endpoint_less{eps}.intersect(b, c, x + R);
            value_type y3 = endpoint_less{eps}.intersect(a, c, x + R);
            value_type y4 = endpoint_less{eps}.intersect(b, a, x + R);
            value_type y5 = endpoint_less{eps}.intersect(c, b, x + R);
            value_type y6 = endpoint_less{eps}.intersect(c, a, x + R);
            asm volatile ("nop");
        }
        return vertices_.insert({{x, y}, R});
    }

    pedge
    add_edge(site const l, site const r, pvertex const v)
    {
        return edges_.insert(std::cend(edges_), {l, r, v, nov});
    }

    pedge
    add_edge(site const l, site const r)
    {
        return add_edge(l, r, nov);
    }

    void
    trunc_edge(edge & e, pvertex const v) const
    {
        assert(v != nov);
        if (e.b == nov) {
            if (e.e == nov) { // orientate if needed:
                point const & l = *e.l;
                point const & r = *e.r;
                point const & c = v->c;
                if (r.x < l.x) {
                    if (c.y < l.y) {
                        e.b = v;
                        return;
                    }
                } else if (l.x < r.x) {
                    if (r.y < c.y) {
                        e.b = v;
                        return;
                    }
                } else {
                    assert(!(r.y < l.y));
                }
                e.e = v;
            } else {
                assert(e.e != v);
                e.b = v;
            }
        } else {
            assert(e.b != v);
            assert(e.e == nov);
            e.e = v;
        }
    }

    void
    disable_event(pvertex const v)
    {
        assert(!events_.empty());
        pevent l = events_.lower_bound(v);
        do {
            pvertex & vv = l->second->second;
            if (vv == v) {
                vv = nov;
            } else {
                assert(vv == nov);
            }
            events_.erase(l++);
        } while ((l != noe) && (l->first == v));
        vertices_.erase(v);
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l);
        static int i = 0;
        std::cerr << ++i << std::endl;
        if (i == 4131) {
            asm volatile ("nop");
        }
        auto const v = make_vertex(*ll.first.l, *ll.first.r, *rr.first.r);
        if (v.first != nov) {
            {
                if (endpoint_less{eps}(*v.first, ll.first)) {
                    endpoint_less{eps}.print_delta(*v.first, ll.first);
                }
                if (endpoint_less{eps}(ll.first, *v.first)) {
                    endpoint_less{eps}.print_delta(*v.first, ll.first);
                }
                if (endpoint_less{eps}(*v.first, rr.first)) {
                    endpoint_less{eps}.print_delta(*v.first, rr.first);
                }
                if (endpoint_less{eps}(rr.first, *v.first)) {
                    endpoint_less{eps}.print_delta(*v.first, rr.first);
                }
                assert(!endpoint_less{eps}(*v.first, ll.first) && !endpoint_less{eps}(ll.first, *v.first));
                assert(!endpoint_less{eps}(*v.first, rr.first) && !endpoint_less{eps}(rr.first, *v.first));
            }
            if ((ll.second != nov) && (rr.second != nov)) {
                value_type const & lx = ll.second->x();
                value_type const & rx = rr.second->x();
                if (lx + eps < rx) {
                    disable_event(rr.second);
                } else if (rx + eps < lx) {
                    disable_event(ll.second);
                }
            }
            if (ll.second != nov) {
                value_type const & lx = v.first->x();
                value_type const & rx = ll.second->x();
                if (lx + eps < rx) {
                    disable_event(ll.second);
                } else {
                    if (rx + eps < lx) {
                        if (v.second) {
                            vertices_.erase(v.first);
                        }
                    } else {
                        assert(!v.second);
                        if (rr.second != v.first) {
                            assert(rr.second == nov);
                            rr.second = v.first;
                            events_.insert({v.first, r});
                        }
                    }
                    return;
                }
            } else if (rr.second != nov) {
                value_type const & lx = v.first->x();
                value_type const & rx = rr.second->x();
                if (lx + eps < rx) {
                    disable_event(rr.second);
                } else {
                    if (rx + eps < lx) {
                        if (v.second) {
                            vertices_.erase(v.first);
                        }
                    } else {
                        assert(!v.second);
                        if (ll.second != v.first) {
                            assert(ll.second == nov);
                            ll.second = v.first;
                            events_.insert({v.first, l});
                        }
                    }
                    return;
                }
            }
            assert(ll.second == nov);
            assert(rr.second == nov);
            ll.second = rr.second = v.first;
            events_.insert(events_.insert({v.first, r}), {v.first, l});
        }
    }

    pendpoint
    insert_endpoint(pendpoint const ep, site const l, site const r, pedge const e)
    {
        return endpoints_.insert(ep, {{l, r, e}, nov});
    }

    void
    finish_endpoints(pendpoint l,
                     pendpoint const r,
                     pvertex const v,
                     site const s)
    {
        assert(1 < std::distance(l, r));
        site const lp = l->first.l;
        site const rp = std::prev(r)->first.r;
        for (pendpoint ep = l; ep != r; ++ep) {
            assert(ep->second == v);
            trunc_edge(*ep->first.e, v);
        }
        events_.erase(v);
        endpoints_.erase(l, r);
        pedge const le = add_edge(lp, s, v);
        pedge const re = add_edge(s, rp, v);
        l = insert_endpoint(r, lp, s, le);
        pendpoint const ep = insert_endpoint(r, s, rp, re);
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != noep) {
            check_event(ep, r);
        }
    }

    void
    make_first_edge(site const l, site const r)
    {
        assert(endpoints_.empty());
        pedge const e = add_edge(l, r);
        pendpoint const le = insert_endpoint(noep, l, r, e);
        if (l->x + eps < r->x)  {
            pendpoint const re = insert_endpoint(noep, r, l, e);
            assert(std::next(le) == re);
        }
    }

    void
    begin_cell(site const s)
    {
        auto lr = endpoints_.equal_range(*s);
        if (lr.first == lr.second) {
            if (lr.first == noep) {
                lr.first = std::prev(noep);
                site const l = lr.first->first.r;
                pedge const e = add_edge(l, s);
                lr.second = insert_endpoint(noep, l, s, e);
                if (l->x + eps < s->x)  {
                    insert_endpoint(noep, s, l, e);
                }
                check_event(lr.first, lr.second);
            } else if (lr.first == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const r = lr.second->first.l;
                pedge const e = add_edge(s, r);
                insert_endpoint(lr.second, r, s, e);
                lr.first = insert_endpoint(lr.second, s, r, e);
                check_event(lr.first, lr.second);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --lr.first;
                site const c = lr.first->first.r;
                assert(c == lr.second->first.l);
                pedge const e = add_edge(c, s);
                pendpoint const l = insert_endpoint(lr.second, c, s, e);
                pendpoint const r = insert_endpoint(lr.second, s, c, e);
                assert(std::next(l) == r);
                check_event(lr.first, l);
                check_event(r, lr.second);
            }
        } else { // many arc collapsing right here, event (equivalent to the current site p) coming on the next step
            finish_endpoints(lr.first, lr.second, lr.first->second, s);
        }
    }

    void
    finish_cells(pevent ev, pvertex const v)
    {
        pendpoint r = ev->second;
        do {
            events_.erase(ev++);
        } while ((ev != noe) && (ev->first == v));
        pendpoint l = r;
        pendpoint const ll = std::begin(endpoints_);
        while (l != ll) {
            if (std::prev(l)->second != v) {
                break;
            }
            --l;
        }
        while (++r != noep) {
            if (r->second != v) {
                break;
            }
        }
        assert(1 < std::distance(l, r));
        site const lc = l->first.l;
        site const rc = std::prev(r)->first.r;
        for (pendpoint ep = l; ep != r; ++ep) {
            assert(ep->second == v);
            trunc_edge(*ep->first.e, v);
        }
        endpoints_.erase(l, r);
        pedge const e = add_edge(lc, rc, v);
        l = insert_endpoint(r, lc, rc, e);
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != noep) {
            check_event(l, r);
        }
    }

    bool
    prior(vertex const & l, point const & r) const
    {
        value_type const & x = l.x() + eps;
        value_type const & y = l.y() + eps;
        return std::tie(x, y) < std::tie(r.x, r.y);
    }

public :

    void
    operator () (site l, site const r)
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
                pevent const e = std::begin(events_);
                pvertex const v = e->first;
                if (!prior(*v, *l)) {
                    break;
                }
                finish_cells(e, v);
            }
            begin_cell(l);
        }
        while (!events_.empty()) {
            pevent const e = std::begin(events_);
            finish_cells(e, e->first);
        }
        endpoints_.clear();
    }

};
