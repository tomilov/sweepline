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
    // (b == inf) means (b == (-infty, infty)); (e == inf) means (e == (+infty, infty))
    // in all other cases (b->c < e->c)
    struct edge
    {

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::size_type;

    vertices vertices_; // size <= 2 * n − 5
    const pvertex inf = std::numeric_limits< pvertex >::max();
    edges edges_; // size <= 3 * n − 6

private :

    template< typename type >
    struct range { type l, r; };

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

        bool operator () (const point & l, const point & r) const
        {
            return operator () (l.x, l.y, r.x, r.y);
        }

        bool operator () (const point & l, const point & r, const point & p, const bool b) const
        {
            const auto f = [this] (const value_type & l, const value_type & r, const bool b) -> bool
            {
                if (b) {
                    return operator () (l, r);
                } else {
                    return operator () (r, l);
                }
            };
            const auto g = [&] (const bool b) -> bool
            {
                value_type dx = (l.x + r.x) / value_type(2) + (p.y - (l.y + r.y) / value_type(2)) * (r.y - l.y) / (l.x - r.x);
                const value_type dxy = p.x - dx;
                dx -= r.x;
                const value_type dy = r.y - p.y;
                return f(dy * dy + dx * dx, dxy * dxy, b);
            };
            if (operator () (l.x, r.x)) {
                if (operator () (p.y, r.y)) {
                    return g(b);
                } else {
                    return b;
                }
            } else if (operator () (r.x, l.x)) {
                if (operator () (l.y, p.y)) {
                    return g(!b);
                } else {
                    return !b;
                }
            } else {
                assert(operator () (l.y, r.y));
                return f((l.y + r.y) / value_type(2), p.y, b);
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
        value_type D = ca.y * cb.x - ca.x * cb.y;
        if (!(less_.eps < D)) { // if not CW
            return {}; // interesting, that probability of this branch is exactly = 0.4 for points in general positions
        }
        D += D;
        const value_type A = ca.x * ca.x + ca.y * ca.y;
        const value_type B = cb.x * cb.x + cb.y * cb.y;
        const value_type x = (B * ca.y - A * cb.y) / D;
        const value_type y = (cb.x * A - ca.x * B) / D;
        using std::sqrt; // std::sqrt is required by the IEEE standard be exact (error < 0.5 ulp)
        return {{{c.x + x, c.y + y}, sqrt(x * x + y * y)}};
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
            pevent ev = events_.search(vertex_, less_);
            const value_type & x = event_x(vertex_);
            const auto deselect_event = [&] (const pevent _ev) -> bool
            {
                if (_ev != nev) {
                    if (_ev != ev) {
                        const value_type & xx = event_x(_ev->k);
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
        pendpoint l, r;
        std::tie(l, r) = endpoints_.equal_range(*s);
        if (l == r) {
            if (l == nep) { // append to the rightmost endpoint
                --l;
                r = add_cell(l->k.r, s);
            } else if (l == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                const site c = r->k.l;
                const pedge e = add_edge(s, c, inf);
                const pendpoint ll = insert_endpoint(r, c, s, e);
                l = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --l;
                const site c = l->k.r;
                assert(c == r->k.l);
                const pedge e = add_edge(c, s, inf);
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
                assert(less_(s->x, event_x(endpoint_.v->k)));
                disable_event(endpoint_.v);
            }
            auto vertex_ = make_vertex(*s, *endpoint_.k.l, *endpoint_.k.r);
            assert(!!vertex_);
            assert(events_.find(*vertex_) == nev);
            assert(!less_(s->x, s->y, event_x(*vertex_), vertex_->c.y)); // vertex and site are equivalent
            assert(!less_(event_x(*vertex_), vertex_->c.y, s->x, s->y)); // vertex and site are equivalent
            const pvertex v = vertices_.size();
            vertices_.push_back(std::move(*vertex_));
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
        const point & point_ = *l;
        while (!events_.empty()) {
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
