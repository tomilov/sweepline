/* Fortune's algorithm implementation
 *
 * Copyright (c) 2016, Anatoliy V. Tomilov
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following condition is met:
 * Redistributions of source code must retain the above copyright notice, this condition and the following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

#include "rb_tree.hpp"

#include <type_traits>
#include <utility>
#include <tuple>
#include <functional>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <deque>
#include <list>
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

    static_assert(std::is_base_of< std::forward_iterator_tag, typename std::iterator_traits< site >::iterator_category >::value,
                  "multipass guarantee required");

    static_assert(std::is_same< decltype(std::declval< point >().x), decltype(std::declval< point >().y) >::value,
                  "point format error");

    sweepline(value_type const && _eps) = delete;

    explicit
    sweepline(value_type const & _eps)
        : less_{_eps}
    {
        assert(!(_eps < value_type(0)));
    }

    struct vertex // circumscribed circle
    {

        point c; // circumcenter
        value_type R; // circumradius

        value_type x() const { return c.x + R; }
        value_type const & y() const { return c.y; }

    };

    using vertices = std::list< vertex >;
    using pvertex = typename vertices::iterator;

    struct edge // ((l, r), (b, e)) is CW
    {

        site l;
        site r;

        pvertex b;
        pvertex e;

        void flip() // !(e->c < b->c) && ((b != nv) || (e == nv))
        {
            using std::swap;
            swap(l, r);
            swap(b, e);
        }

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::size_type;

    // Voronoi diagram:
    // NOTE: logically diagram is neither copyable nor moveable due to past the end iterator of std::list is not preserved during these operations
    // {
    vertices vertices_;
    pvertex const nv = std::end(vertices_); // infty
    edges edges_;
    // }

private :

    static
    value_type
    angle(point const & l, point const & r)
    {
        using std::atan2;
        return atan2(r.x - l.x, r.y - l.y);
    }

    struct endpoint
    {

        site const l;
        site const r;

        pedge const e;

        value_type angle() const { return sweepline::angle(*l, *r); }

    };

    struct less
    {

        value_type const & eps;

        bool operator () (value_type const & l,
                          value_type const & r) const
        {
            return l + eps < r;
        }

        bool operator () (value_type const & lx, value_type const & ly,
                          value_type const & rx, value_type const & ry) const
        {
            if (operator () (lx, rx)) {
                return true;
            } else if (operator () (rx, lx)) {
                return false;
            } else {
                return operator () (ly, ry);
            }
        }

        bool operator () (point const & l, point const & r) const
        {
            return operator () (l.x, l.y, r.x, r.y);
        }

        bool operator () (vertex const & l, vertex const & r) const
        {
            return operator () (l.x(), l.y(), r.x(), r.y());
        }

        // During sweepline motion arcs shrinks and growz, but relative y-order of endpoints remains the same.
        // Endpoints removed strictly before a violation of this invariant to prevent its occurrence.

        // The body of this function is unreachable.
        // It still gives correct results, but only when used against std::map from libc++, but doing so is UB.
        //[[noreturn]]
        bool operator () (endpoint const & l, endpoint const & r) const
        {
            //__builtin_unreachable();
            //throw this;
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
            return operator () (std::max(*l.l, *l.r, *this), std::max(*r.l, *r.r, *this));
        }

        value_type
        intersect(point const & l, point const & r,
                  value_type const & directrix) const
        {
            if (operator () (r.x, l.x)) {
                if (!operator () (l.x, directrix)) {
                    assert(!operator () (directrix, l.x));
                    return l.y;
                }
            } else if (operator () (l.x, r.x)) {
                if (!operator () (r.x, directrix)) {
                    assert(!operator () (directrix, r.x));
                    return r.y;
                }
            } else { // a ~= 0
                assert(!operator () (directrix, l.x));
                assert(!operator () (directrix, r.x));
                return (l.y + r.y) / value_type(2);
            }
            value_type const a = l.x - r.x;
            value_type const ld = l.x - directrix;
            value_type const rd = r.x - directrix;
            value_type const b = r.y * ld - l.y * rd; // -b
            value_type const c = r.y * r.y * ld - l.y * l.y * rd - ld * rd * a;
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

        bool operator () (point const & l, endpoint const & r) const
        {
            return operator () (l.y, intersect(r, l.x));
        }

        bool operator () (endpoint const & l, point const & r) const
        {
            return operator () (intersect(l, r.x), r.y);
        }

    } const less_;

    struct pevent;

    using endpoints = rb_tree::map< endpoint, pevent, less >;
    using pendpoint = typename endpoints::iterator;

    using rays = std::list< pendpoint >;
    using pray = typename rays::iterator;

    template< typename type >
    struct range { type l, r; };

    using bundle = range< pray const >;

    using events = rb_tree::map< vertex, bundle const, less >;

    using event = typename events::iterator;
    struct pevent : event { pevent(event const it) : event{it} { ; } };

    endpoints endpoints_{less_};
    pendpoint const nep = std::end(endpoints_);

    rays rays_;
    pray const nray = std::end(rays_);
    pray rev = nray; // revocation boundary

    events events_{less_};
    pevent const nev = std::end(events_);

    std::experimental::optional< vertex >
    make_vertex(point const & a,
                point const & b,
                point const & c) const
    {
        point const ca = {a.x - c.x, a.y - c.y};
        point const cb = {b.x - c.x, b.y - c.y};
        value_type alpha = ca.x * cb.y - ca.y * cb.x;
        if (!(less_.eps < -alpha)) {
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
        assert(less_.eps * less_.eps < beta + x * x + y * y);
        using std::sqrt; // std::sqrt is required by the IEEE standard be exact (error < 0.5 ulp)
        return {{{x, y}, sqrt(beta + x * x + y * y)}};
    }

    void add_ray(pray const rr, pendpoint const l)
    {
        assert(rr != nray);
        if (nray == rev) {
            rays_.insert(rr, l);
        } else {
            *rev = l;
            rays_.splice(rr, rays_, rev++);
        }
    }

    bundle
    add_bundle(pendpoint const l, pendpoint const r)
    {
        if (rev == nray) {
            return {rays_.insert(nray, l), rays_.insert(nray, r)};
        } else {
            pray const ll = rev;
            *ll = l;
            if (++rev == nray) {
                return {ll, rays_.insert(nray, r)};
            } else {
                *rev = r;
                return {ll, rev++};
            }
        }
    }

    void remove_bundle(bundle const & b)
    {
        rays_.splice(nray, rays_, b.l, std::next(b.r));
        if (rev == nray) {
            rev = b.l;
        }
    }

    void disable_event(pevent const ev)
    {
        assert(ev != nev);
        bundle const & b = ev->v;
        assert(b.l != b.r);
        assert(nray != b.r);
        remove_bundle(b);
        assert(std::next(b.r) == nray);
        for (auto l = b.l; l != nray; ++l) {
            pendpoint const ep = *l;
            assert(ep->v == ev);
            ep->v = nev;
        }
        events_.erase(ev);
    }

    void check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.k.r == rr.k.l);
        if (auto v = make_vertex(*ll.k.l, *ll.k.r, *rr.k.r)) {
            vertex & vertex_ = *v;
            pevent ev = events_.find(vertex_);
            value_type const & x = vertex_.x();
            auto const deselect_event = [&] (pevent const _ev) -> bool
            {
                if (_ev != nev) {
                    if (_ev != ev) {
                        value_type const & xx = _ev->k.x();
                        if (less_(xx, x)) {
                            return true;
                        }
                        assert(less_(x, xx)); // not equiv
                        disable_event(_ev);
                    }
                }
                return false;
            };
            if (deselect_event(ll.v) || deselect_event(rr.v)) {
                if (ev != nev) {
                    disable_event(ev);
                }
            } else {
                if (ev == nev) {
                    assert(ll.v == nev);
                    assert(rr.v == nev);
                    bool inserted = false;
                    std::tie(ev, inserted) = events_.insert({std::move(vertex_), add_bundle(l, r)});
                    assert(inserted);
                    ll.v = rr.v = ev;
                } else {
                    bundle const & b = ev->v;
                    auto const set_event = [&] (pevent & _ev, pendpoint const ep)
                    {
                        if (_ev == nev) {
                            _ev = ev;
                            add_ray(b.r, ep);
                        } else {
                            assert(_ev == ev);
                            assert(std::find(b.l, std::next(b.r), ep) != std::next(b.r));
                        }
                    };
                    set_event(ll.v, l);
                    set_event(rr.v, r);
                }
            }
        }
    }

    pedge
    add_edge(site const l, site const r, pvertex const v)
    {
        assert(l != r);
        pedge const e = edges_.size();
        edges_.push_back({l, r, v, nv});
        return e;
    }

    void trunc_edge(pedge const e, pvertex const v)
    {
        edge & edge_ = edges_[e];
        assert(v != nv);
        if (edge_.b == nv) {
            if (edge_.e == nv) { // orientate:
                point const & l = *edge_.l;
                point const & r = *edge_.r;
                point const & c = v->c;
                if (r.x < l.x) {
                    if (c.y < l.y) {
                        edge_.b = v;
                        return;
                    }
                } else if (l.x < r.x) {
                    if (r.y < c.y) {
                        edge_.b = v;
                        return;
                    }
                } else {
                    assert(!(r.y < l.y));
                }
                edge_.e = v;
                edge_.flip();
                return;
            } else {
                assert(edge_.e != v);
                edge_.b = v;
            }
        } else {
            assert(edge_.b != v);
            assert(edge_.e == nv);
            edge_.e = v;
        }
        point const & l = edge_.b->c;
        point const & r = edge_.e->c;
        if (std::tie(r.x, r.y) < std::tie(l.x, l.y)) {
            edge_.flip();
        }
    }

    pendpoint
    insert_endpoint(pendpoint const ep,
                    site const l, site const r,
                    pedge const e)
    {
        return endpoints_.force_insert(ep, {{l, r, e}, nev});
    }

    pendpoint
    add_cell(site const c, site const s)
    {
        pedge const e = add_edge(c, s, nv);
        pendpoint const r = insert_endpoint(nep, c, s, e);
        if (less_(c->x, s->x))  {
            pendpoint const rr = insert_endpoint(nep, s, c, e);
            assert(std::next(r) == rr);
        }
        return r;
    }

    void begin_cell(site const s)
    {
        assert(!endpoints_.empty());
        pendpoint l = endpoints_.lower_bound(*s);
        pendpoint r = l;
        while ((r != nep) && !less_(*s, r->k)) {
            assert(l->v == r->v); // if fires, then there is problem with precision
            ++r;
        }
        if (l == r) {
            if (l == nep) { // append to the rightmost endpoint
                --l;
                r = add_cell(l->k.r, s);
            } else if (l == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const c = r->k.l;
                pedge const e = add_edge(s, c, nv);
                pendpoint const ll = insert_endpoint(r, c, s, e);
                l = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --l;
                site const c = l->k.r;
                assert(c == r->k.l);
                pedge const e = add_edge(c, s, nv);
                pendpoint const ll = insert_endpoint(r, c, s, e);
                pendpoint const rr = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == rr);
                check_event(l, ll);
                check_event(rr, r);
                return;
            }
            check_event(l, r);
        } else {
            assert(std::next(l) == r); // if fires, then there is problem with precision
            auto const & endpoint_ = *l;
            if (endpoint_.v != nev) {
                assert(less_(s->x, endpoint_.v->k.x()));
                disable_event(endpoint_.v);
            }
            auto vertex_ = make_vertex(*s, *endpoint_.k.l, *endpoint_.k.r);
            assert(!!vertex_);
            assert(events_.find(*vertex_) == nev);
            pvertex const v = vertices_.insert(nv, std::move(*vertex_));
            trunc_edge(endpoint_.k.e, v);
            pedge const le = add_edge(endpoint_.k.l, s, v);
            pedge const re = add_edge(s, endpoint_.k.r, v);
            pendpoint const ep = insert_endpoint(r, s, endpoint_.k.r, re);
            assert(std::next(ep) == r);
            endpoints_.erase(std::exchange(l, insert_endpoint(ep, endpoint_.k.l, s, le)));
            assert(std::next(l) == ep);
            if (l != std::begin(endpoints_)) {
                check_event(std::prev(l), l);
            }
            if (r != nep) {
                check_event(ep, r);
            }
        }
    }

    range< pendpoint >
    endpoint_range(pray l, pray r)
    {
        assert(r != nray);
        assert(l != r);
        if (std::next(l) == r) {
            assert(std::next(*l) == *r);
        } else {
            auto const angle_less = [] (pendpoint const ll, pendpoint const rr) -> bool
            {
                return ll->k.angle() < rr->k.angle();
            };
#if 0
            assert(std::next(r) == nray);
            rays crays_;
            crays_.splice(std::cend(crays_), rays_, l, nray);
            crays_.sort(angle_less);
            l = std::begin(crays_);
            rays_.splice(nray, std::move(crays_));
            r = std::prev(nray);
#else
            std::tie(l, r) = std::minmax_element(l, std::next(r), angle_less);
#endif
        }
        return {*l, *r};
    }

    bool check_endpoint_range(pevent const ev, pendpoint l, pendpoint const r) const
    {
        if (r == nep) {
            return false;
        }
        if (l == r) {
            return false;
        }
        site s = l->k.r;
        do {
            ++l;
            if (std::exchange(s, l->k.r) != l->k.l) {
                return false;
            }
            if (l->v != ev) {
                return false;
            }
        } while (l != r);
        return true;
    }

    void finish_cells(pevent const ev,
                      vertex const & _vertex,
                      bundle const & b,
                      site const l, site const r)
    {
        remove_bundle(b);
        auto lr = endpoint_range(b.l, b.r);
        // All the edges from (*(lr.l ... lr.r))->k.e can be stored near the associate vertex if needed
        assert(check_endpoint_range(ev, lr.l, lr.r));
        pvertex const v = vertices_.insert(nv, _vertex);
        events_.erase(ev);
        site const ll = lr.l->k.l;
        site const rr = lr.r->k.r;
        ++lr.r;
        do {
            trunc_edge(lr.l->k.e, v);
            endpoints_.erase(lr.l++);
        } while (lr.l != lr.r);
        if (l == r) {
            lr.l = insert_endpoint(lr.r, ll, rr, add_edge(ll, rr, v));
            if (lr.l != std::begin(endpoints_)) {
                check_event(std::prev(lr.l), lr.l);
            }
            if (lr.r != nep) {
                check_event(lr.l, lr.r);
            }
        } else {
            pendpoint const ep = insert_endpoint(lr.r, l, rr, add_edge(l, rr, v));
            lr.l = insert_endpoint(ep, ll, l, add_edge(ll, l, v));
            assert(std::next(lr.l) == ep);
            if (lr.l != std::begin(endpoints_)) {
                check_event(std::prev(lr.l), lr.l);
            }
            if (lr.r != nep) {
                check_event(ep, lr.r);
            }
        }
    }

    bool check_last_endpoints() const
    {
        for (auto const & ep : endpoints_) {
            edge const & e = edges_[ep.k.e];
            if ((e.b != nv) && (e.e != nv)) {
                return false;
            }
            if (ep.v != nev) {
                return false;
            }
        }
        return true;
    }

    bool process_events(site const l, site const r)
    {
        point const & point_ = *l;
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            auto const & event_ = *ev;
            value_type const & x = event_.k.x();
            if (less_(point_.x, x)) {
                break;
            } else if (!less_(x, point_.x)) {
                value_type const & y = event_.k.y();
                if (less_(point_.y, y)) {
                    break;
                } else if (!less_(y, point_.y)) {
                    finish_cells(ev, event_.k, event_.v, l, r);
                    return false;
                }
            }
            finish_cells(ev, event_.k, event_.v, l, l);
        }
        return true;
    }

public :

    template< typename iterator >
    void operator () (iterator l, iterator const r)
    {
        assert(std::is_sorted(l, r, less_));
        assert(endpoints_.empty());
        assert(vertices_.empty());
        assert(edges_.empty());
        if (l == r) {
            return;
        }
        iterator const ll = std::next(l);
        if (ll == r) {
            return;
        }
        add_cell(l, ll);
        l = ll;
        while (++l != r) {
            if (process_events(l, r)) {
                begin_cell(l);
            }
        }
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            auto const & event_ = *ev;
            finish_cells(ev, event_.k, event_.v, r, r);
        }
        //assert(std::is_sorted(std::begin(vertices_), nv, less_)); // almost true
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
