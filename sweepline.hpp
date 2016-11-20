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
    bool // imprecise comparator
    less(value_type const & l,
         value_type const & _eps,
         value_type const & r)
    {
        return l + _eps < r;
    }

    static
    bool // lexicographically compare with tolerance
    less(value_type const & lx, value_type const & ly,
         value_type const & _eps,
         value_type const & rx, value_type const & ry)
    {
        if (less(lx, _eps, rx)) {
            return true;
        } else if (less(rx, _eps, lx)) {
            return false;
        } else {
            if (less(ly, _eps, ry)) {
                return true;
            } else {
                return false;
            }
        }
    }

    struct vertex // circumscribed circle
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
            return less(l.x, l.y, eps_, r.x, r.y);
        }

        bool operator () (site const l, site const r) const
        {
            return operator () (*l, *r);
        }

    };

    using vertices = std::list< vertex >;
    using pvertex = typename vertices::iterator;

    struct edge // ((l, r), (b, e)) is CW
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::iterator;

    // Voronoi diagram:
    // NOTE: logically diagram is neither copyable nor moveable due to past the end iterator of std::list is not preserved during this operations
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
    intersect(point const & focus,
              value_type const & y,
              value_type const & directrix)
    {
        assert(!(directrix/* + eps*/ < focus.x));
        value_type d = focus.x - directrix;
        return (y * (y - (focus.y + focus.y)) + (focus.x * focus.x + focus.y * focus.y - directrix * directrix)) / (d + d);
    }

    struct endpoint_less
    {

        value_type const & eps_;

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
            point_less const point_less_{eps_};
            return point_less_(std::max(*l.l, *l.r, point_less_), std::max(*r.l, *r.r, point_less_));
        }

        using is_transparent = void;

        value_type
        intersect(point const & l,
                  point const & r,
                  value_type const & directrix) const
        {
            if (less(r.x, eps_, l.x)) {
                if (!less(l.x, eps_, directrix)) {
                    assert(!less(directrix, eps_, l.x));
                    return l.y;
                }
            } else if (less(l.x, eps_, r.x)) {
                if (!less(r.x, eps_, directrix)) {
                    assert(!less(directrix, eps_, r.x));
                    return r.y;
                }
            } else { // a ~= 0
                assert(!less(directrix, eps_, l.x));
                assert(!less(directrix, eps_, r.x));
                return (l.y + r.y) / value_type(2);
            }
            value_type const a = l.x - r.x;
            value_type const ld = l.x - directrix;
            value_type const rd = r.x - directrix;
            value_type const b = r.y * ld - l.y * rd; // -b
            value_type const c = (rd * r.x + r.y * r.y) * ld - (ld * l.x + l.y * l.y) * rd;
            value_type const D = b * b - a * c;
            assert(!(D < value_type(0)));
            using std::sqrt;
            return (b + sqrt(D)) / a;
        }

        value_type
        intersect(endpoint const & ep, value_type const & directrix) const
        {
            return intersect(*ep.l, *ep.r, directrix);
        }

        bool operator () (vertex const & l, endpoint const & r) const
        {
            return less(l.y(), eps_, intersect(r, l.x()));
        }

        bool operator () (endpoint const & l, vertex const & r) const
        {
            return less(intersect(l, r.x()), eps_, r.y());
        }

        bool operator () (point const & l, endpoint const & r) const
        {
            return less(l.y, eps_, intersect(r, l.x));
        }

        bool operator () (endpoint const & l, point const & r) const
        {
            return less(intersect(l, r.x), eps_, r.y);
        }

    };

    struct event;
    using endpoints = std::map< endpoint, event, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    struct event_less
    {

        value_type const & eps_;

        bool operator () (vertex const & l, vertex const & r) const
        {
            return less(l.x(), l.y(), eps_, r.x(), r.y());
        }

    };

    using rays = std::list< pendpoint >;
    using pray = typename rays::iterator;

    using bundle = std::pair< pray const, pray const >;

    using events = std::map< vertex, bundle const, event_less >;
    using pevent = typename events::iterator;

    struct event
    {

        pevent ev;

        event(pevent const _ev) : ev(_ev) { ; }

        operator pevent () const { return ev; }
        operator pevent & () { return ev; }

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

    std::experimental::optional< vertex >
    make_vertex(point const & a,
                point const & b,
                point const & c) const
    {
        point const ca = {a.x - c.x, a.y - c.y};
        point const cb = {b.x - c.x, b.y - c.y};
        value_type alpha = ca.x * cb.y - ca.y * cb.x;
        if (!(eps < -alpha)) {
            return {};
        }
        value_type const A = a.x * a.x + a.y * a.y;
        value_type const B = b.x * b.x + b.y * b.y;
        value_type const C = c.x * c.x + c.y * c.y;
        value_type const CA = A - C;
        value_type const CB = B - C;
        value_type x = CA * cb.y - CB * ca.y;
        value_type y = ca.x * CB - cb.x * CA;
        value_type beta = a.x * (b.y * C - c.y * B) + b.x * (c.y * A - a.y * C) + c.x * (a.y * B - b.y * A);
        beta /= alpha;
        alpha += alpha;
        x /= alpha;
        y /= alpha;
        assert(eps * eps < beta + x * x + y * y);
        using std::sqrt; // std::sqrt is required by the IEEE standard be exact (error < 0.5 ulp)
        return {{{x, y}, sqrt(beta + x * x + y * y)}};
    }

    void
    trunc_edge(edge & e, pvertex const v) const
    {
        assert(v != nov);
        if (e.b == nov) {
            if (e.e == nov) { // orientate:
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
    add_ray(pray const rr, pendpoint const l)
    {
        assert(rr != noray);
        if (noray == rev) {
            rays_.insert(rr, l);
        } else {
            *rev = l;
            rays_.splice(rr, rays_, rev++);
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
        assert(b.first != b.second);
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
            pevent ev = events_.find(*v);
            value_type const & x = v->x();
            auto const deselect_event = [&] (pevent const _ev) -> bool
            {
                if (_ev != noev) {
                    if (_ev != ev) {
                        value_type const & xx = _ev->first.x();
                        if (less(xx, eps, x)) {
                            return true;
                        }
                        assert(less(x, eps, xx));
                        disable_event(_ev);
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
                    bundle const & b = ev->second;
                    auto const set_event = [&] (pevent & _ev, pendpoint const ep)
                    {
                        if (_ev == noev) {
                            _ev = ev;
                            add_ray(b.second, ep);
                        } else {
                            assert(_ev == ev);
                            assert(std::find(b.first, std::next(b.second), ep) != std::next(b.second));
                        }
                    };
                    set_event(ll.second, l);
                    set_event(rr.second, r);
                }
            }
        }
    }

    pendpoint
    insert_endpoint(pendpoint const ep,
                    site const l, site const r,
                    pedge const e)
    {
        return endpoints_.insert(ep, {{l, r, e}, noev});
    }

    void
    make_first_edge(site const l, site const r)
    {
        assert(endpoints_.empty());
        pedge const e = add_edge(l, r, nov);
        pendpoint const ll = insert_endpoint(noep, l, r, e);
        if (less(l->x, eps, r->x))  {
            pendpoint const rr = insert_endpoint(noep, r, l, e);
            assert(std::next(ll) == rr);
        }
    }

    void
    finish_endpoints(pendpoint l,
                     pendpoint const r,
                     site const s)
    {
        endpoint const & endpoint_ = l->first;
        auto vertex_ = make_vertex(*s, *endpoint_.l, *endpoint_.r);
        assert(!!vertex_);
        assert(events_.find(*vertex_) == noev);
        pvertex const v = vertices_.insert(nov, std::move(*vertex_));
        trunc_edge(*endpoint_.e, v);
        pedge const ll = add_edge(endpoint_.l, s, v);
        pedge const rr = add_edge(s, endpoint_.r, v);
        pendpoint const ep = insert_endpoint(r, s, endpoint_.r, rr);
        assert(std::next(ep) == r);
        endpoints_.erase(std::exchange(l, insert_endpoint(ep, endpoint_.l, s, ll)));
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
        pendpoint l = endpoints_.lower_bound(*s);
        pendpoint r = l;
        while (r != noep) {
            if (endpoint_less{eps}(*s, r->first)) {
                break;
            }
            assert(l->second.ev == r->second.ev);
            ++r;
        }
        if (l == r) {
            if (l == noep) { // append to the rightmost endpoint
                --l;
                site const c = l->first.r;
                pedge const e = add_edge(c, s, nov);
                r = insert_endpoint(noep, c, s, e);
                if (less(c->x, eps, s->x))  {
                    pendpoint const rr = insert_endpoint(noep, s, c, e);
                    assert(std::next(r) == rr);
                }
            } else if (l == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const c = r->first.l;
                pedge const e = add_edge(s, c, nov);
                pendpoint const ll = insert_endpoint(r, c, s, e);
                l = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --l;
                site const c = l->first.r;
                assert(c == r->first.l);
                pedge const e = add_edge(c, s, nov);
                pendpoint const ll = insert_endpoint(r, c, s, e);
                pendpoint const rr = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == rr);
                check_event(l, ll);
                check_event(rr, r);
                return;
            }
            check_event(l, r);
        } else {
            assert(std::next(l) == r); // problems with precision detected
            if (l->second != noev) {
                assert(less(s->x, eps, l->second.ev->first.x()));
                disable_event(l->second);
            }
            finish_endpoints(l, r, s);
        }
    }

    static
    value_type
    angle(point const & l, point const & r)
    {
        using std::atan2;
        return atan2(r.x - l.x, r.y - l.y);
    }

    std::pair< pendpoint, pendpoint >
    endpoint_range(pray const l, pray const r) const
    {
        assert(r != noray);
        assert(0 < std::distance(l, r));
        if (std::next(l) == r) {
            assert(std::next(*l) == *r);
            return {*l, *r};
        } else { // [l; r] can be sorted (to be stored with associate vertex) right here if needed in some application
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

    bool
    check_endpoint_range(pevent const ev, pendpoint l, pendpoint const r) const
    {
        if (r == noep) {
            return false;
        }
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
                 bundle const & b,
                 std::experimental::optional< site const > const _site = {})
    {
        auto lr = endpoint_range(b.first, b.second);
        assert(check_endpoint_range(ev, lr.first, lr.second));
        pvertex const v = vertices_.insert(nov, _vertex);
        remove_event(ev, b);
        site const l = lr.first->first.l;
        site const r = lr.second->first.r;
        ++lr.second;
        do {
            trunc_edge(*lr.first->first.e, v);
            endpoints_.erase(lr.first++);
        } while (lr.first != lr.second);
        if (_site) {
            site const s = *_site;
            pedge const le = add_edge(l, s, v);
            pedge const re = add_edge(s, r, v);
            pendpoint const ep = insert_endpoint(lr.second, s, r, re);
            lr.first = insert_endpoint(ep, l, s, le);
            assert(std::next(lr.first) == ep);
            if (lr.first != std::begin(endpoints_)) {
                check_event(std::prev(lr.first), lr.first);
            }
            if (lr.second != noep) {
                check_event(ep, lr.second);
            }
        } else {
            pedge const e = add_edge(l, r, v);
            lr.first = insert_endpoint(lr.second, l, r, e);
            if (lr.first != std::begin(endpoints_)) {
                check_event(std::prev(lr.first), lr.first);
            }
            if (lr.second != noep) {
                check_event(lr.first, lr.second);
            }
        }
    }

    bool
    prior(vertex const & l, point const & r) const
    {
        return less(l.x(), l.y(), eps, r.x, r.y);
    }

    bool
    prior(point const & l, vertex const & r) const
    {
        return less(l.x, l.y, eps, r.x(), r.y());
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
        assert(std::is_sorted(l, r, point_less{eps}));
        assert(vertices_.empty());
        assert(edges_.empty());
        if (l == r) {
            return;
        }
        if (std::next(l) == r) {
            return;
        }
        make_first_edge(l, ++l);
        while (++l != r) {
            bool continue_ = false;
            while (!events_.empty()) {
                pevent const ev = std::begin(events_);
                auto & event_ = *ev;
                if (!prior(event_.first, *l)) {
                    if (!prior(*l, event_.first)) {
                        finish_cells(ev, event_.first, event_.second, l);
                        continue_ = true;
                    }
                    break;
                }
                finish_cells(ev, event_.first, event_.second);
            }
            if (!continue_) {
                begin_cell(l);
            }
        }
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            auto & event_ = *ev;
            finish_cells(ev, event_.first, event_.second);
        }
        assert(rev == std::begin(rays_));
        assert(check_last_endpoints());
        endpoints_.clear();
    }

    void clear()
    {
        assert(rev == std::begin(rays_));
        assert(endpoints_.empty());
        assert(events_.empty());
        vertices_.clear();
        edges_.clear();
    }

};
