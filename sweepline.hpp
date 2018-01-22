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

    };

    using vertices = std::deque< vertex >;
    using pvertex = typename vertices::size_type;

    // ((l, r), (b, e)) is CW
    // (b == inf) means (b == (-infty, *)); (e == inf) means (e == (+infty, *))
    // in all other cases (b->c < e->c)
    struct edge
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::size_type;

    vertices vertices_; // 0 <= size <= 2 * n − 5
    const pvertex inf = std::numeric_limits< pvertex >::max();
    edges edges_; // n - 1 <= size <= 3 * n − 6

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

    static
    value_type event_x(const vertex & v)
    {
        return v.c.x + v.R;
    }

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
            return operator () (event_x(l), l.c.y, event_x(r), r.c.y);
        }

        bool operator () (const point & l, const point & r, const point & p, const bool right) const
        {
            const auto swap_cmp = [this] (const value_type & ll, const value_type & rr, const bool b) -> bool
            {
                if (b) {
                    return operator () (ll, rr);
                } else {
                    return operator () (rr, ll);
                }
            };
            const auto sqr_dist = [&] (const bool b) -> bool
            {
                const point c = {(l.x + r.x) / value_type(2), (l.y + r.y) / value_type(2)};
                value_type dx = c.x + (p.y - c.y) * (r.y - l.y) / (l.x - r.x);
                const value_type dxy = p.x - dx;
                dx -= r.x;
                const value_type dy = r.y - p.y;
                return swap_cmp(dy * dy + dx * dx, dxy * dxy, b);
            };
            if (operator () (l.x, r.x)) {
                if (operator () (p.y, r.y)) {
                    return sqr_dist(right);
                } else {
                    return right;
                }
            } else if (operator () (r.x, l.x)) {
                if (operator () (l.y, p.y)) {
                    return sqr_dist(!right);
                } else {
                    return !right;
                }
            } else {
                assert(operator () (l.y, r.y));
                return swap_cmp((l.y + r.y) / value_type(2), p.y, right);
            }
        }

        bool operator () (const point & l, const endpoint & r) const
        {
            return operator () (*r.l, *r.r, l, false);
        }

        bool operator () (const endpoint & l, const point & r) const
        {
            return operator () (*l.l, *l.r, r, true);
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

    vertex make_vertex(const point & a,
                       const point & b,
                       const point & c) const
    {
        const point ca = {a.x - c.x, a.y - c.y};
        const point cb = {b.x - c.x, b.y - c.y};
        vertex vertex_{{}, ca.y * cb.x - ca.x * cb.y};
        if (less_.eps < vertex_.R) { // if CW
            // interesting, that probability of this branch tends to 0.6 for points in general positions
            vertex_.R += vertex_.R;
            const value_type A = ca.x * ca.x + ca.y * ca.y;
            const value_type B = cb.x * cb.x + cb.y * cb.y;
            vertex_.c.x = (B * ca.y - A * cb.y) / vertex_.R;
            vertex_.c.y = (cb.x * A - ca.x * B) / vertex_.R;
            using std::sqrt; // std::sqrt is required by the IEEE standard be exact (error < 0.5 ulp)
            vertex_.R = sqrt(vertex_.c.x * vertex_.c.x + vertex_.c.y * vertex_.c.y);
            vertex_.c.x += c.x;
            vertex_.c.y += c.y;
        }
        return vertex_;
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

    pevent disable_event(const pevent ev)
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
        return events_.erase(ev);
    }

    void check_event(const pendpoint l, const pendpoint r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.k.r == rr.k.l);
        vertex vertex_ = make_vertex(*ll.k.l, *ll.k.r, *rr.k.r);
        if (less_.eps < vertex_.R) {
            const value_type & x = event_x(vertex_);
            const auto le = events_.find(vertex_);
            const auto deselect_event = [&] (const pevent ev) -> bool
            {
                if (ev != nev) {
                    if (ev != le) {
                        const value_type & xx = event_x(ev->k);
                        if (less_(xx, x)) {
                            return true;
                        }
                        assert(less_(x, xx)); // not equiv
                        disable_event(ev);
                    }
                }
                return false;
            };
            if (deselect_event(ll.v) || deselect_event(rr.v)) {
                if (le != nev) {
                    disable_event(le);
                }
            } else {
                if (le == nev) {
                    assert(ll.v == nev);
                    assert(rr.v == nev);
                    const auto ev = events_.insert({std::move(vertex_), add_bundle(l, r)});
                    assert(ev.v);
                    ll.v = rr.v = ev.k;
                } else {
                    const bundle & b = le->v;
                    const auto set_event = [&] (pevent & ev, const pendpoint ep)
                    {
                        if (ev == nev) {
                            ev = le;
                            add_ray(b.r, ep);
                        } else {
                            assert(ev == le);
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
        const point & ll = *l;
        const point & rr = *r;
        const pedge e = edges_.size();
        if ((ll.y < rr.y) || (!(rr.y < ll.y) && (rr.x < ll.x))) {
            edges_.push_back({l, r, v, inf});
        } else {
            edges_.push_back({r, l, inf, v});
        }
        return e;
    }

    void truncate_edge(const pedge e, const pvertex v)
    {
        assert(v != inf);
        edge & edge_ = edges_[e];
        if (edge_.b != inf) {
            assert(edge_.b != v);
            assert(edge_.e == inf);
            edge_.e = v;
        } else if (edge_.e != inf) {
            assert(edge_.e != v);
            edge_.b = v;
        } else {
            const point & l = *edge_.l;
            const point & r = *edge_.r;
            const point & c = vertices_[v].c;
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
        // workaround for floating point math
        point & ve = vertices_[edge_.e].c;
        point & vb = vertices_[edge_.b].c;
        if (ve.x < vb.x) {
            ve.x = vb.x = (vb.x + ve.x) / value_type(2);
        }
        assert(vb < ve);
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
        const pedge e = add_edge(c, s, inf);
        const pendpoint r = insert_endpoint(nep, c, s, e);
        if (less_(c->x, s->x))  {
            const pendpoint rr = insert_endpoint(nep, s, c, e);
            assert(std::next(r) == rr);
        }
        return r;
    }

    void begin_cell(const site s)
    {
        assert(!endpoints_.empty());
        auto lr = endpoints_.equal_range(*s);
        if (lr.l == lr.r) {
            if (lr.l == nep) { // append to the rightmost endpoint
                --lr.l;
                lr.r = add_cell(lr.l->k.r, s);
            } else if (lr.l == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                const site c = lr.r->k.l;
                const pedge e = add_edge(s, c, inf);
                const pendpoint ll = insert_endpoint(lr.r, c, s, e);
                lr.l = insert_endpoint(lr.r, s, c, e);
                assert(std::next(ll) == lr.l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --lr.l;
                const site c = lr.l->k.r;
                assert(c == lr.r->k.l);
                const pedge e = add_edge(c, s, inf);
                const pendpoint ll = insert_endpoint(lr.r, c, s, e);
                const pendpoint rr = insert_endpoint(lr.r, s, c, e);
                assert(std::next(ll) == rr);
                check_event(lr.l, ll);
                check_event(rr, lr.r);
                return;
            }
            check_event(lr.l, lr.r);
        } else {
            assert(std::next(lr.l) == lr.r); // if fires, then there is problem with precision
            const auto & endpoint_ = *lr.l;
            if (endpoint_.v != nev) {
                assert(less_(s->x, event_x(endpoint_.v->k)));
                disable_event(endpoint_.v);
            }
            vertex vertex_ = make_vertex(*s, *endpoint_.k.l, *endpoint_.k.r);
            assert(less_.eps < vertex_.R);
            assert(events_.find(vertex_) == nev);
            assert(!less_(s->x, s->y, event_x(vertex_), vertex_.c.y)); // vertex and site are equivalent
            assert(!less_(event_x(vertex_), vertex_.c.y, s->x, s->y)); // vertex and site are equivalent
            const pvertex v = vertices_.size();
            vertices_.push_back(std::move(vertex_));
            truncate_edge(endpoint_.k.e, v);
            const pedge le = add_edge(endpoint_.k.l, s, v);
            const pedge re = add_edge(s, endpoint_.k.r, v);
            const pendpoint ep = insert_endpoint(lr.r, s, endpoint_.k.r, re);
            assert(std::next(ep) == lr.r);
            endpoints_.erase(std::exchange(lr.l, insert_endpoint(ep, endpoint_.k.l, s, le)));
            assert(std::next(lr.l) == ep);
            if (lr.l != std::begin(endpoints_)) {
                check_event(std::prev(lr.l), lr.l);
            }
            if (lr.r != nep) {
                check_event(ep, lr.r);
            }
        }
    }

    range< pendpoint >
    endpoint_range(pray l, pray r)
    {
        assert(l != nray);
        assert(r != nray);
        assert(std::next(r) == nray);
        assert(l != r);
        if (std::next(l) == r) {
            assert(std::next(*l) == *r);
        } else {
            const auto angle_less = [] (const pendpoint ll, const pendpoint rr) -> bool
            {
                return ll->k.angle() < rr->k.angle();
            };
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
                assert(false);
                return false;
            }
            if (l->v != ev) {
                assert(false);
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
        assert(check_endpoint_range(ev, lr.l, lr.r));
        const pvertex v = vertices_.size();
        vertices_.push_back(_vertex);
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
            if ((e.b != inf) && (e.e != inf)) {
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
        if (!events_.empty()) {
            const point & point_ = *l;
            do {
                const pevent ev = std::begin(events_);
                const auto & event_ = *ev;
                const value_type & x = event_x(event_.k);
                if (less_(point_.x, x)) {
                    break;
                } else if (!less_(x, point_.x)) {
                    const value_type & y = event_.k.c.y;
                    if (less_(point_.y, y)) {
                        break;
                    } else if (!less_(y, point_.y)) {
                        finish_cells(ev, event_.k, event_.v, l, r);
                        return false;
                    }
                }
                finish_cells(ev, event_.k, event_.v, l, l);
            } while (!events_.empty());
        }
        return true;
    }

public :

    template< typename iterator >
    void operator () (iterator l, const iterator r)
    {
        static_assert(std::is_base_of< std::forward_iterator_tag, typename std::iterator_traits< iterator >::iterator_category >::value,
                      "multipass guarantee required");
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
