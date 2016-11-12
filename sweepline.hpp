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
#include <list>
#include <set>
#include <map>
#include <experimental/optional>
#ifdef DEBUG
#include <iostream>
#endif

#include <cassert>
#include <cmath>

template< typename site,
          typename point = typename std::iterator_traits< site >::value_type,
          typename value_type = decltype(std::declval< point >().x) >
struct sweepline
{

    static_assert(std::is_same< decltype(std::declval< point >().x), decltype(std::declval< point >().y) >::value,
                  "point format error");

    value_type const & eps;

    sweepline(value_type const && _eps) = delete;

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    {
        assert(!(eps < value_type(0)));
    }

    static
    bool // imprecise comparator: lexicographically compare with tolerance
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

        bool operator () (site const l, site const r) const
        {
            return operator () (*l, *r);
        }

    };

    using vertices = std::list< vertex >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is CCW
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::iterator;

    // Voronoi diagram:
    // NOTE: logically diagram is neither copyable nor moveable due to past the end iterator is not preserved during this operations
    // {
    vertices vertices_;
    pvertex const nov = std::end(vertices_);
    edges edges_;
    // }

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
        assert(!(directrix/* + eps*/ < p.x));
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
            /*if (&l == &r) {
                return false;
            }*/
            // All the next code up to "}" is (striclty saying) undefined behaviour,
            // but it works for all modern standard libraries and (I sure)
            // more stricter then current local (insertion by hint) guarantees
            // for comp() should be the part of the Standard
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
            point_less point_less_{eps_};
            return point_less_(std::max(*l.l, *l.r, point_less_), std::max(*r.l, *r.r, point_less_));
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

    struct wevent;
    using endpoints = std::map< endpoint, wevent, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    struct event_less
    {

        value_type const & eps_;

        bool operator () (vertex const & l, vertex const & r) const
        {
            return less(eps_, l.x(), l.y(), r.x(), r.y());
        }

    };

    using rays = std::list< pendpoint >;
    using pray = typename rays::iterator;

    using bundle = std::pair< pray, pray >;

    using events = std::map< vertex, bundle, event_less >;
    using pevent = typename events::iterator;

    struct wevent
    {

        wevent(pevent const _ev)
            : ev(_ev)
        { ; }

        operator pevent () const
        {
            return ev;
        }

        operator pevent & ()
        {
            return ev;
        }

        pevent ev;

    };

    endpoints endpoints_{endpoint_less{eps}};
    pendpoint const noep = std::end(endpoints_);

    rays rays_;
    pray const noray = std::end(rays_);
    pray rev = noray; // revocation boundary

    events events_{event_less{eps}};
    pevent const noev = std::end(events_);

    pedge
    add_edge(site const l, site const r, pvertex const v)
    {
        return edges_.insert(std::cend(edges_), {l, r, v, nov});
    }

    void
    make_first_edge(site const l, site const r)
    {
        assert(endpoints_.empty());
        pedge const e = add_edge(l, r, nov);
        pendpoint const le = insert_endpoint(noep, l, r, e);
        if (l->x + eps < r->x)  {
            pendpoint const re = insert_endpoint(noep, r, l, e);
            assert(std::next(le) == re);
        }
    }

    value_type
    circumradius(value_type const a,
                 value_type const b,
                 value_type const c) const
    {
        value_type const V = (a + b - c) * (a + c - b) * (b + c - a);
        assert(eps * eps * eps < V + V); // triangle inequality
        using std::sqrt;
        return (a * b * c) / sqrt(V * (a + b + c));
    }

    std::experimental::optional< vertex >
    make_vertex(point const & a,
                point const & b,
                point const & c)
    {
        value_type const A = a.x * a.x + a.y * a.y;
        value_type const B = b.x * b.x + b.y * b.y;
        value_type const C = c.x * c.x + c.y * c.y;
        point const ca = {a.x - c.x, a.y - c.y};
        point const cb = {b.x - c.x, b.y - c.y};
        value_type x = (A - C) * cb.y - (B - C) * ca.y;
        value_type y = ca.x * (B - C) - cb.x * (A - C);
        value_type alpha = ca.x * cb.y - ca.y * cb.x;
        if (!(eps * eps < -alpha)) {
            return {};
        }
        value_type beta = a.x * (b.y * C - c.y * B) - b.x * (a.y * C - c.y * A) + c.x * (a.y * B - b.y * A);
        beta /= alpha;
        alpha += alpha;
        x /= alpha;
        y /= alpha;
        using std::sqrt; // std::sqrt is required by the IEEE standard be exact (error < 0.5 ulp)
        value_type R = sqrt(beta + x * x + y * y);
        return {{{x, y}, R}};
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
    add_ray(pray const r, pendpoint const ep)
    {
        assert(r != noray);
        if (noray == rev) {
            rays_.insert(r, ep);
        } else {
            *rev = ep;
            rays_.splice(r, rays_, rev++);
        }
    }

    bundle
    add_bundle(pendpoint const l, pendpoint const r)
    {
        if (rev == noray) {
            return {rays_.insert(noray, l), rays_.insert(noray, r)};
        } else {
            pray const ll = rev;
            *ll = l;
            if (++rev == noray) {
                return {ll, rays_.insert(noray, r)};
            } else {
                *rev = r;
                return {ll, rev++};
            }
        }
    }

    void
    remove_event(pevent const ev, bundle const & b)
    {
        rays_.splice(noray, rays_, b.first, std::next(b.second));
        if (rev == noray) {
            rev = b.first;
        }
        events_.erase(ev);
    }

    void
    disable_event(pevent const ev)
    {
        assert(ev != noev);
        bundle const & b = ev->second;
        pray const r = std::next(b.second);
        for (auto l = b.first; l != r; ++l) {
            pendpoint const ep = *l;
            assert(ep->second == ev);
            ep->second = noev;
        }
        remove_event(ev, b);
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l);
        if (auto v = make_vertex(*ll.first.l, *ll.first.r, *rr.first.r)) {
            auto ev = events_.find(*v);
            value_type const & x = v->x();
            auto const deselect_event = [&] (pevent const pev) -> bool
            {
                if (pev != noev) {
                    if (pev != ev) {
                        value_type const & xx = pev->first.x();
                        if (xx + eps < x) {
                            return true;
                        }
                        assert(x + eps < xx);
                        disable_event(pev);
                    }
                }
                return false;
            };
            if (deselect_event(ll.second) || deselect_event(rr.second)) {
                if (ev != noev) {
                    disable_event(ev);
                }
            } else {
                if (ev == noev) {
                    assert(ll.second == noev);
                    assert(rr.second == noev);
                    bool inserted = false;
                    std::tie(ev, inserted) = events_.insert({std::move(*v), add_bundle(l, r)});
                    assert(inserted);
                    ll.second = rr.second = ev;
                } else {
                    bundle & b = ev->second;
                    auto const set_event = [&] (pevent & pev, pendpoint const lr)
                    {
                        if (pev == noev) {
                            pev = ev;
                            add_ray(b.second, lr);
                        } else {
                            assert(pev == ev);
                            assert(std::find(b.first, std::next(b.second), lr) != std::next(b.second));
                        }
                    };
                    set_event(rr.second, r);
                    set_event(ll.second, l);
                }
            }
        }
    }

    pendpoint
    insert_endpoint(pendpoint const ep,
                    site const l, site const r,
                    pedge const e)
    {
        return endpoints_.insert(ep, {{l, r, e}, {noev}});
    }

    void
    finish_endpoints(pendpoint l,
                     pendpoint const r,
                     pevent ev,
                     site const s)
    {
        site const lf = l->first.l;
        site const rf = std::prev(r)->first.r;
        auto const create_vertex = [&] () -> vertex
        {
            if (ev != noev) {
                if (std::next(l) == r) {
                    assert(s->x + eps < ev->first.x());
                    disable_event(ev);
                    ev = noev;
                } else {
                    assert(!(s->x + eps < ev->first.x()));
                    return ev->first;
                }
            }
            assert(std::next(l) == r);
            auto v = make_vertex(*s, *lf, *rf);
            assert(!!v);
            return std::move(*v);
        };
        pvertex const v = vertices_.insert(nov, create_vertex());
        do {
            assert(l->second == ev);
            trunc_edge(*l->first.e, v);
            endpoints_.erase(l++);
        } while (l != r);
        if (ev != noev) {
            remove_event(ev, ev->second);
        }
        pedge const le = add_edge(lf, s, v);
        pedge const re = add_edge(s, rf, v);
        pendpoint const ep = insert_endpoint(r, s, rf, re);
        l = insert_endpoint(ep, lf, s, le);
        assert(std::next(l) == ep);
        if (l != std::begin(endpoints_)) {
            check_event(std::prev(l), l);
        }
        if (r != noep) {
            check_event(ep, r);
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
                pedge const e = add_edge(l, s, nov);
                lr.second = insert_endpoint(noep, l, s, e);
                if (l->x + eps < s->x)  {
                    pendpoint const ep = insert_endpoint(noep, s, l, e);
                    assert(std::next(lr.second) == ep);
                }
            } else if (lr.first == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const r = lr.second->first.l;
                pedge const e = add_edge(s, r, nov);
                insert_endpoint(lr.second, r, s, e);
                lr.first = insert_endpoint(lr.second, s, r, e);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --lr.first;
                site const c = lr.first->first.r;
                assert(c == lr.second->first.l);
                pedge const e = add_edge(c, s, nov);
                pendpoint const l = insert_endpoint(lr.second, c, s, e);
                pendpoint const r = insert_endpoint(lr.second, s, c, e);
                assert(std::next(l) == r);
                check_event(lr.first, l);
                check_event(r, lr.second);
                return;
            }
            check_event(lr.first, lr.second);
        } else { // one endpoint or many arc collapsing right here, event (equivalent to the current site p) coming on the next step
            { // workaround for libc++ bug https://llvm.org/bugs/show_bug.cgi?id=30959
                auto const ll = std::begin(endpoints_);
                while (lr.first != ll) {
                    --lr.first;
                    if (endpoint_less{eps}(lr.first->first, *s)) {
                        ++lr.first;
                        break;
                    }
                }
                if (lr.second != noep) {
                    while (!endpoint_less{eps}(*s, lr.second->first)) {
                        if (++lr.second == noep) {
                            break;
                        }
                    }
                }
            }
            finish_endpoints(lr.first, lr.second, lr.first->second, s);
        }
    }

    static
    value_type
    angle(point const & l, point const & r)
    {
        using std::atan2;
        return atan2(r.x - l.x, r.y - l.y);
    }

    static
    std::pair< pendpoint, pendpoint >
    boundaries(pray const l, pray const r)
    {
        assert(0 < std::distance(l, r));
        if (std::next(l) == r) {
            assert(std::next(*l) == *r);
            return {*l, *r};
        } else { // [l; r] can be sorted right here if needed in some application
            auto const angle_less = [&] (pendpoint const ll, pendpoint const rr) -> bool
            {
                endpoint const & lll = ll->first;
                endpoint const & rrr = rr->first;
                return angle(*lll.l, *lll.r) < angle(*rrr.l, *rrr.r);
            };
            auto const lr = std::minmax_element(l, std::next(r), angle_less);
            return {*lr.first, *lr.second};
        }
    }

    static
    bool
    check_endpoint_range(pevent const ev, pendpoint l, pendpoint const r)
    {
        if (l == r) {
            return false;
        }
        site s = l->first.r;
        do {
            ++l;
            if (std::exchange(s, l->first.r) != l->first.l) {
                return false;
            }
            if (l->second != ev) {
                return false;
            }
        } while (l != r);
        return true;
    }

    void
    finish_cells(pevent const ev,
                 vertex const & _vertex,
                 bundle const & b)
    {
        auto lr = boundaries(b.first, b.second);
        assert(check_endpoint_range(ev, lr.first, lr.second));
        pvertex const v = vertices_.insert(nov, _vertex);
        remove_event(ev, b);
        site const lc = lr.first->first.l;
        site const rc = lr.second->first.r;
        ++lr.second;
        do {
            trunc_edge(*lr.first->first.e, v);
            endpoints_.erase(lr.first++);
        } while (lr.first != lr.second);
        pedge const e = add_edge(lc, rc, v);
        lr.first = insert_endpoint(lr.second, lc, rc, e);
        if (lr.first != std::begin(endpoints_)) {
            check_event(std::prev(lr.first), lr.first);
        }
        if (lr.second != noep) {
            check_event(lr.first, lr.second);
        }
    }

    bool
    prior(vertex const & l, point const & r) const
    {
        return less(eps, l.x(), l.y(), r.x, r.y);
    }

    bool
    check_last_endpoints() const
    {
        for (auto const & ep : endpoints_) {
            edge const & e = *ep.first.e;
            if ((e.b != nov) && (e.e != nov)) {
                return false;
            }
            if (ep.second != noev) {
                return false;
            }
        }
        return true;
    }

public :

    template< typename iterator >
    void
    operator () (iterator l, iterator const r)
    {
        assert(vertices_.empty());
        assert(edges_.empty());
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
                auto & event_ = *ev;
                if (!prior(event_.first, *l)) {
                    break;
                }
                finish_cells(ev, event_.first, event_.second);
            }
            begin_cell(l);
        }
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            auto & event_ = *ev;
            finish_cells(ev, event_.first, event_.second);
        }
        assert(check_last_endpoints());
        assert(rev == std::begin(rays_));
        endpoints_.clear();
        rays_.clear();
    }

    void clear()
    {
        vertices_.clear();
        edges_.clear();
    }

};
