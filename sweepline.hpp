#pragma once

#include <type_traits>
#include <utility>
#include <tuple>
#include <functional>
#include <iterator>
#include <limits>
#include <algorithm>
#include <numeric>
#include <deque>
#include <set>
#include <map>

#include <cassert>
#include <cmath>

template< typename site,
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

    static
    bool
    less(value_type const & _eps,
         value_type const & lx, value_type const & ly,
         value_type const & rx, value_type const & ry)
    {
        if (lx + _eps < rx) {
            return true;
        } else if (rx + _eps < lx) {
            return false;
        } else {
            if (ly + _eps < ry) {
                return true;
            } else {
                return false;
            }
        }
    }

    struct vertex // denote circumscribed circle
    {

        point c; // circumcenter
        value_type R; // circumradius

        value_type x() const { return c.x + R; }
        value_type const & y() const { return c.y; }

    };

    struct point_less
    {

        value_type const & eps_;

        bool operator () (point const & l, point const & r) const
        {
            return less(eps_, l.x, l.y, r.x, r.y);
        }

        bool operator () (vertex const & l, vertex const & r) const
        {
            return operator () (l.c, r.c);
        }

    };

    using vertices = std::set< vertex, point_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::iterator;

    vertices vertices_{point_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;

private :

    struct endpoint
    {

        site l, r;
        pedge e;

    };

    static
    value_type
    intersect(point const & p,
              value_type const & y,
              value_type const & directrix)
    {
        //assert(!(directrix + eps < p.x));
        value_type d = p.x - directrix;
        return (y * (y - (p.y + p.y)) + (p.x * p.x + p.y * p.y - directrix * directrix)) / (d + d);
    }

    struct endpoint_less
    {

        value_type const & eps_;

        value_type
        intersect(point const & l,
                  point const & r,
                  value_type directrix) const
        {
            {
                bool const rdegenerated = !(r.x + eps_ < directrix);
                if (!(l.x + eps_ < directrix)) {
                    assert(!(directrix + eps_ < l.x));
                    if (rdegenerated) {
                        assert(!(directrix + eps_ < r.x));
                        assert((l.y + eps_ < r.y) || (r.y + eps_ < l.y)); // l != r
                        return (l.y + r.y) / value_type(2);
                    } else {
                        return l.y;
                    }
                } else if (rdegenerated) {
                    assert(!(directrix + eps_ < r.x));
                    return r.y;
                }
            }
            value_type ld = l.x - directrix;
            value_type rd = r.x - directrix;
            value_type b = r.y / rd - l.y / ld; // -b
            ld += ld;
            rd += rd;
            directrix *= directrix;
            value_type lc = (l.x * l.x + l.y * l.y - directrix) / ld;
            value_type rc = (r.x * r.x + r.y * r.y - directrix) / rd;
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
        intersect(endpoint const & ep, value_type const & directrix) const
        {
            return intersect(*ep.l, *ep.r, directrix);
        }

        // during sweepline motion arcs shrinks and growz, but relative y-position of endpoints remains the same
        // endpoints removed strictly before violation of this invariant to prevent its occurrence
        bool operator () (endpoint const & l, endpoint const & r) const
        {
            // It is undefined behaviour to use following four "if"s,
            // but it works for all modern standard libraries and I sure
            // local (insertion by hint) guarantees for comp() should be the part of the Standard
            if (l.l == r.l) {
                return true;
            }
            if (r.r == l.r) {
                return true;
            }
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            assert(false);
            if (&l == &r) {
                return false;
            }
            point const & ll = *l.l;
            point const & lr = *l.r;
            point const & rl = *r.l;
            point const & rr = *r.r;
            value_type const & directrix = std::max(std::max(ll.x, lr.x), std::max(rl.x, rr.x));
            return intersect(ll, lr, directrix) + eps_ < intersect(rl, rr, directrix);
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
            return less(eps_, l.x(), l.y(), r.x(), r.y());
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

    value_type
    circumradius(value_type a,
                 value_type b,
                 value_type c) const
    {
        value_type V = (a + b - c) * (a + c - b) * (b + c - a);
        assert(eps * eps * eps < V); // triangle inequality
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
        value_type M = A * (a.x + b.x) + B * (a.y + b.y);
        value_type N = C * (b.x + c.x) + D * (b.y + c.y);
        G += G;
        // circumcenter:
        value_type y = (C * M - A * N) / G;
        auto const miss = [&] (point const & l, point const & r) -> bool
        {
            if (l.x < r.x) {
                if (r.y < y) {
                    return true;
                }
            } else {
                if (y < l.y) {
                    return true;
                }
            }
            return false;
        };
        if (miss(a, b) || miss(b, c)) {
            return {nov, false};
        }
        value_type x = (B * N - D * M) / G;
        using std::hypot;
#if 1
        value_type R = circumradius(hypot(A, B), hypot(C, D), hypot(a.x - c.x, a.y - c.y));
#else
        value_type R = hypot(x - b.x, y - b.y);
#endif
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
        assert(v != nov);
        assert(!events_.empty());
        pevent ev = events_.lower_bound(v);
        assert(ev != noe);
        assert(ev->first == v);
        do {
            assert(ev->second->second == v);
            ev->second->second = nov;
            events_.erase(ev++);
        } while ((ev != noe) && (ev->first == v));
        vertices_.erase(v);
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l);
        auto const v = make_vertex(*ll.first.l, *ll.first.r, *rr.first.r);
        if (v.first != nov) {
            assert((events_.find(v.first) == noe) == v.second);
            value_type const & x = v.first->x();
            auto const deselect_event = [&] (auto const & ep) -> bool
            {
                if (ep.second != nov) {
                    if (ep.second != v.first) {
                        value_type const & xx = ep.second->x();
                        if (xx + eps < x) {
                            if (v.second) {
                                vertices_.erase(v.first);
                            } else {
                                disable_event(v.first);
                            }
                            return false;
                        }
                        assert(x + eps < xx);
                        disable_event(ep.second);
                    }
                }
                return true;
            };
            if (!deselect_event(ll)) {
                return;
            }
            if (!deselect_event(rr)) {
                return;
            }
            auto const set_event = [&] (auto & ep, pendpoint const lr)
            {
                if (ep.second == nov) {
                    ep.second = v.first;
                    events_.insert({v.first, lr});
                } else {
                    assert(ep.second == v.first);
                }
            };
            set_event(ll, l);
            set_event(rr, r);
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
                     pvertex v,
                     site const s)
    {
        site const lp = l->first.l;
        site const rp = std::prev(r)->first.r;
        if (v == nov) {
            assert(std::next(l) == r);
            bool inserted = false;
            std::tie(v, inserted) = make_vertex(*s, *lp, *rp);
            assert(inserted);
            l->second = v;
        } else {
            assert(events_.find(v) != noe);
            events_.erase(v);
        }
        do {
            assert(l->second == v);
            trunc_edge(*l->first.e, v);
            endpoints_.erase(l++);
        } while (l != r);
        pedge const le = add_edge(lp, s, v);
        pedge const re = add_edge(s, rp, v);
        if (r == noep) {
            l = insert_endpoint(r, lp, s, le);
            insert_endpoint(r, s, rp, re);
        } else {
            pendpoint const ep = insert_endpoint(r, s, rp, re);
            l = insert_endpoint(ep, lp, s, le);
            check_event(ep, r);
        }
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
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
            } else if (lr.first == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const r = lr.second->first.l;
                pedge const e = add_edge(s, r);
                insert_endpoint(lr.second, r, s, e);
                lr.first = insert_endpoint(lr.second, s, r, e);
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
                return;
            }
            check_event(lr.first, lr.second);
        } else { // many arc collapsing right here, event (equivalent to the current site p) coming on the next step
            finish_endpoints(lr.first, lr.second, lr.first->second, s);
        }
    }

    void
    finish_cells(pevent ev, pvertex const v)
    {
        pendpoint r = ev->second;
        do {
            assert(ev->second->second == v);
            events_.erase(ev++);
        } while ((ev != noe) && (ev->first == v));
        pendpoint l = r;
        pendpoint const ll = std::begin(endpoints_);
        while (l != ll) {
            if ((--l)->second != v) {
                ++l;
                break;
            }
        }
        while (++r != noep) {
            if (r->second != v) {
                break;
            }
        }
        assert(1 < std::distance(l, r));
        site const lc = l->first.l;
        site const rc = std::prev(r)->first.r;
        do {
            assert(l->second == v);
            trunc_edge(*l->first.e, v);
            endpoints_.erase(l++);
        } while (l != r);
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
        return less(eps, l.x(), l.y(), r.x, r.y);
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
                pevent const ev = std::begin(events_);
                pvertex const v = ev->first;
                if (!prior(*v, *l)) {
                    break;
                }
                finish_cells(ev, v);
            }
            begin_cell(l);
        }
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            finish_cells(ev, ev->first);
        }
        endpoints_.clear();
    }

};
