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

    static_assert(std::is_same< decltype(std::declval< point >() < std::declval< point >()), bool >::value,
                  "points cannot be ordered");

    sweepline(const value_type && _eps) = delete;

    explicit
    sweepline(const value_type & _eps)
        : less_{_eps}
    {
        assert(!(_eps < value_type(0)));
    }

    struct vertex // circumscribed circle
    {

        point c; // circumcenter
        value_type R; // circumradius

        value_type x() const { return c.x + R; }
        const value_type & y() const { return c.y; }

    };

    using vertices = std::list< vertex >;
    using pvertex = typename vertices::iterator;

    // ((l, r), (b, e)) is CW
    // if (b == nv), then (b == (-infty, infty)), if (e == nv), then (e == (+infty, infty))
    // in all other cases (b->c < e->c)
    struct edge
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::size_type;

    // Voronoi diagram:
    // NOTE: logically diagram is neither copyable nor moveable due to past the end iterator of std::list is not preserved during these operations
    // {
    vertices vertices_;
    const pvertex nv = std::end(vertices_); // infty
    edges edges_;
    // }

private :

    struct endpoint
    {

        const site l;
        const site r;

        const pedge e;

        value_type angle() const
        {
            const point & ll = *l;
            const point & rr = *r;
            using std::atan2;
            return atan2(rr.x - ll.x, rr.y - ll.y);
        }

    };

    struct less
    {

        const value_type & eps;

        bool operator () (const value_type & l,
                          const value_type & r) const
        {
            return l + eps < r;
        }

        bool operator () (const value_type & lx, const value_type & ly,
                          const value_type & rx, const value_type & ry) const
        {
            if (operator () (lx, rx)) {
                return true;
            } else if (operator () (rx, lx)) {
                return false;
            } else {
                return operator () (ly, ry);
            }
        }

        bool operator () (const vertex & l, const vertex & r) const
        {
            return operator () (l.x(), l.y(), r.x(), r.y());
        }

        bool operator () (const endpoint & l, const endpoint & r) const = delete;

        value_type intersect(const point & l, const point & r,
                             const value_type & directrix) const
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
            const value_type a = l.x - r.x;
            const value_type ld = l.x - directrix;
            const value_type rd = r.x - directrix;
            const value_type b = r.y * ld - l.y * rd; // -b
            const value_type c = r.y * r.y * ld - l.y * l.y * rd - ld * rd * a;
            const value_type D = b * b - a * c;
            assert(!(D < value_type(0)));
            using std::sqrt;
            return (b + sqrt(D)) / a;
        }

        value_type intersect(const endpoint & ep, const value_type & directrix) const
        {
            return intersect(*ep.l, *ep.r, directrix);
        }

        bool operator () (const point & l, const endpoint & r) const
        {
            return operator () (l.y, intersect(r, l.x));
        }

        bool operator () (const endpoint & l, const point & r) const
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

    using bundle = range< const pray >;

    using events = rb_tree::map< vertex, bundle const, less >;

    using pevent_base = typename events::iterator;
    struct pevent : pevent_base { pevent(const pevent_base it) : pevent_base{it} { ; } };

    endpoints endpoints_{less_};
    const pendpoint nep = std::end(endpoints_);

    rays rays_;
    const pray nray = std::end(rays_);
    pray rev = nray; // revocation boundary

    events events_{less_};
    const pevent nev = std::end(events_);

    std::experimental::optional< vertex >
    make_vertex(const point & a,
                const point & b,
                const point & c) const
    {
        const point ca = {a.x - c.x, a.y - c.y};
        const point cb = {b.x - c.x, b.y - c.y};
        value_type alpha = ca.x * cb.y - ca.y * cb.x;
        if (!(less_.eps < -alpha)) {
            return {};
        }
        const value_type A = a.x * a.x + a.y * a.y;
        const value_type B = b.x * b.x + b.y * b.y;
        const value_type C = c.x * c.x + c.y * c.y;
        const value_type CA = A - C;
        const value_type CB = B - C;
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

    void add_ray(const pray rr, const pendpoint l)
    {
        assert(rr != nray);
        if (nray == rev) {
            rays_.insert(rr, l);
        } else {
            *rev = l;
            rays_.splice(rr, rays_, rev++);
        }
    }

    bundle add_bundle(const pendpoint l, const pendpoint r)
    {
        if (rev == nray) {
            return {rays_.insert(nray, l), rays_.insert(nray, r)};
        } else {
            const pray ll = rev;
            *ll = l;
            if (++rev == nray) {
                return {ll, rays_.insert(nray, r)};
            } else {
                *rev = r;
                return {ll, rev++};
            }
        }
    }

    void remove_bundle(const bundle & b)
    {
        rays_.splice(nray, rays_, b.l, std::next(b.r));
        if (rev == nray) {
            rev = b.l;
        }
    }

    void disable_event(const pevent ev)
    {
        assert(ev != nev);
        const bundle & b = ev->v;
        assert(b.l != b.r);
        assert(nray != b.r);
        remove_bundle(b);
        assert(std::next(b.r) == nray);
        for (auto l = b.l; l != nray; ++l) {
            const pendpoint ep = *l;
            assert(ep->v == ev);
            ep->v = nev;
        }
        events_.erase(ev);
    }

    void check_event(const pendpoint l, const pendpoint r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.k.r == rr.k.l);
        if (auto v = make_vertex(*ll.k.l, *ll.k.r, *rr.k.r)) {
            vertex & vertex_ = *v;
            pevent ev = events_.find(vertex_);
            const value_type & x = vertex_.x();
            const auto deselect_event = [&] (const pevent _ev) -> bool
            {
                if (_ev != nev) {
                    if (_ev != ev) {
                        const value_type & xx = _ev->k.x();
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
                    const bundle & b = ev->v;
                    const auto set_event = [&] (pevent & _ev, const pendpoint ep)
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

    pedge add_edge(const site l, const site r, const pvertex v)
    {
        assert(l != r);
        const pedge e = edges_.size();
        const point & ll = *l;
        const point & rr = *r;
        if ((ll.y < rr.y) || (!(rr.y < ll.y) && (rr.x < ll.x))) {
            edges_.push_back({l, r, v, nv});
        } else {
            edges_.push_back({r, l, nv, v});
        }
        return e;
    }

    void truncate_edge(const pedge e, const pvertex v)
    {
        assert(v != nv);
        edge & edge_ = edges_[e];
        if (edge_.b != nv) {
            assert(edge_.b != v);
            assert(edge_.e == nv);
            edge_.e = v;
        } else if (edge_.e != nv) {
            assert(edge_.e != v);
            edge_.b = v;
        } else {
            const point & l = *edge_.l;
            const point & r = *edge_.r;
            const point & c = v->c;
            assert(!(r.y < l.y));
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
                assert(l.y < c.y);
                assert(c.y < r.y);
            }
            edge_.e = v;
            return;
        }
        if (edge_.e->c < edge_.b->c) {
            std::swap(edge_.b, edge_.e);
            std::swap(edge_.l, edge_.r);
        }
    }

    pendpoint insert_endpoint(const pendpoint ep,
                              const site l, const site r,
                              const pedge e)
    {
        return endpoints_.force_insert(ep, {{l, r, e}, nev});
    }

    pendpoint add_cell(const site c, const site s)
    {
        assert(*c < *s);
        const pedge e = add_edge(c, s, nv);
        const pendpoint r = insert_endpoint(nep, c, s, e);
        if (c->x < s->x)  {
            const pendpoint rr = insert_endpoint(nep, s, c, e);
            assert(std::next(r) == rr);
        }
        return r;
    }

    void begin_cell(const site s)
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
                const site c = r->k.l;
                const pedge e = add_edge(s, c, nv);
                const pendpoint ll = insert_endpoint(r, c, s, e);
                l = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --l;
                const site c = l->k.r;
                assert(c == r->k.l);
                const pedge e = add_edge(c, s, nv);
                const pendpoint ll = insert_endpoint(r, c, s, e);
                const pendpoint rr = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == rr);
                check_event(l, ll);
                check_event(rr, r);
                return;
            }
            check_event(l, r);
        } else {
            assert(std::next(l) == r); // if fires, then there is problem with precision
            const auto & endpoint_ = *l;
            if (endpoint_.v != nev) {
                assert(less_(s->x, endpoint_.v->k.x())); // ?
                disable_event(endpoint_.v);
            }
            auto vertex_ = make_vertex(*s, *endpoint_.k.l, *endpoint_.k.r);
            assert(!!vertex_);
            assert(events_.find(*vertex_) == nev);
            const pvertex v = vertices_.insert(nv, std::move(*vertex_));
            truncate_edge(endpoint_.k.e, v);
            const pedge le = add_edge(endpoint_.k.l, s, v);
            const pedge re = add_edge(s, endpoint_.k.r, v);
            const pendpoint ep = insert_endpoint(r, s, endpoint_.k.r, re);
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
            const auto angle_less = [] (const pendpoint ll, const pendpoint rr) -> bool
            {
                return ll->k.angle() < rr->k.angle();
            };
            assert(std::next(r) == nray);
#if 0
            rays crays_;
            crays_.splice(std::cend(crays_), rays_, l, nray);
            crays_.sort(angle_less);
            l = std::begin(crays_);
            rays_.splice(nray, std::move(crays_));
            r = std::prev(nray);
#else
            std::tie(l, r) = std::minmax_element(l, nray, angle_less);
#endif
        }
        return {*l, *r};
    }

    bool check_endpoint_range(const pevent ev, pendpoint l, const pendpoint r) const
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

    void finish_cells(const pevent ev,
                      const vertex & _vertex,
                      const bundle & b,
                      const site l, const site r)
    {
        remove_bundle(b);
        auto lr = endpoint_range(b.l, b.r);
        // All the edges from (*(lr.l ... lr.r))->k.e can be stored near the associate vertex if needed
        assert(check_endpoint_range(ev, lr.l, lr.r));
        const pvertex v = vertices_.insert(nv, _vertex);
        events_.erase(ev);
        const site ll = lr.l->k.l;
        const site rr = lr.r->k.r;
        ++lr.r;
        do {
            truncate_edge(lr.l->k.e, v);
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
            const pendpoint ep = insert_endpoint(lr.r, l, rr, add_edge(l, rr, v));
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
        for (const auto & ep : endpoints_) {
            const edge & e = edges_[ep.k.e];
            if ((e.b != nv) && (e.e != nv)) {
                return false;
            }
            if (ep.v != nev) {
                return false;
            }
        }
        return true;
    }

    bool process_events(const site l, const site r)
    {
        const point & point_ = *l;
        while (!events_.empty()) {
            const pevent ev = std::begin(events_);
            const auto & event_ = *ev;
            const value_type & x = event_.k.x();
            if (less_(point_.x, x)) {
                break;
            } else if (!less_(x, point_.x)) {
                const value_type & y = event_.k.y();
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
    void operator () (iterator l, const iterator r)
    {
        assert(std::is_sorted(l, r));
        assert(endpoints_.empty());
        assert(vertices_.empty());
        assert(edges_.empty());
        if (l == r) {
            return;
        }
        const iterator ll = l;
        if (++l == r) {
            return;
        }
        add_cell(ll, l);
        while (++l != r) {
            if (process_events(l, r)) {
                begin_cell(l);
            }
        }
        while (!events_.empty()) {
            const pevent ev = std::begin(events_);
            const auto & event_ = *ev;
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
