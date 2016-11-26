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

    static_assert(std::is_base_of< std::forward_iterator_tag, typename std::iterator_traits< site >::iterator_category >::value,
                  "multipass guarantee required");

    static_assert(std::is_same< decltype(std::declval< point >().x), decltype(std::declval< point >().y) >::value,
                  "point format error");

    sweepline(sweepline const &) = delete;
    sweepline & operator = (sweepline const &) = delete;

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

        site l, r;
        pvertex b, e;

    };

    using edges = std::deque< edge >;
    using pedge = typename edges::iterator;

    // Voronoi diagram:
    // NOTE: logically diagram is neither copyable nor moveable due to past the end iterator of std::list is not preserved during these operations
    // {
    vertices vertices_;
    pvertex const nv = std::end(vertices_); // inf
    edges edges_;
    // }

private :

    struct endpoint
    {

        site l, r;
        pedge e;

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

        // during sweepline motion arcs shrinks and growz, but relative y-position of endpoints remains the same
        // endpoints removed strictly before violation of this invariant to prevent its occurrence
        bool operator () (endpoint const & l, endpoint const & r) const // UB
        {
            /*if (&l == &r) {
                return false;
            }*/
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

        using is_transparent = void;

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

    struct event;
    using endpoints = std::map< endpoint, event, less >;
    using pendpoint = typename endpoints::iterator;

    using rays = std::list< pendpoint >;
    using pray = typename rays::iterator;

    using bundle = std::pair< pray const, pray const >;

    using events = std::map< vertex, bundle const, less >;
    using pevent = typename events::iterator;

    struct event
    {

        pevent ev;

        event(pevent const _ev) : ev(_ev) { ; }

        operator pevent () const { return ev; }
        operator pevent & () { return ev; }

    };

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

    void remove_event(pevent const ev, bundle const & b)
    {
        rays_.splice(nray, rays_, b.first, std::next(b.second));
        if (rev == nray) {
            rev = b.first;
        }
        events_.erase(ev);
    }

    void disable_event(pevent const ev)
    {
        assert(ev != nev);
        bundle const & b = ev->second;
        assert(b.first != b.second);
        assert(nray != b.second);
        pray const r = std::next(b.second);
        for (auto l = b.first; l != r; ++l) {
            pendpoint const ep = *l;
            assert(ep->second == ev);
            ep->second = nev;
        }
        remove_event(ev, b);
    }

    void check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l);
        if (auto v = make_vertex(*ll.first.l, *ll.first.r, *rr.first.r)) {
            vertex & vertex_ = *v;
            pevent ev = events_.find(vertex_);
            value_type const & x = vertex_.x();
            auto const deselect_event = [&] (pevent const _ev) -> bool
            {
                if (_ev != nev) {
                    if (_ev != ev) {
                        value_type const & xx = _ev->first.x();
                        if (less_(xx, x)) {
                            return true;
                        }
                        assert(less_(x, xx)); // not equiv
                        disable_event(_ev);
                    }
                }
                return false;
            };
            if (deselect_event(ll.second) || deselect_event(rr.second)) {
                if (ev != nev) {
                    disable_event(ev);
                }
            } else {
                if (ev == nev) {
                    assert(ll.second == nev);
                    assert(rr.second == nev);
                    bool inserted = false;
                    std::tie(ev, inserted) = events_.insert({std::move(vertex_), add_bundle(l, r)});
                    assert(inserted);
                    ll.second = rr.second = ev;
                } else {
                    bundle const & b = ev->second;
                    auto const set_event = [&] (pevent & _ev, pendpoint const ep)
                    {
                        if (_ev == nev) {
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

    pedge
    add_edge(site const l, site const r, pvertex const v)
    {
        return edges_.insert(std::cend(edges_), {l, r, v, nv});
    }

    void trunc_edge(edge & e, pvertex const v) const
    {
        assert(v != nv);
        if (e.b == nv) {
            if (e.e == nv) { // orientate:
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
            assert(e.e == nv);
            e.e = v;
        }
    }

    pendpoint
    insert_endpoint(pendpoint const ep,
                    site const l, site const r,
                    pedge const e)
    {
        return endpoints_.insert(ep, {{l, r, e}, nev});
    }

    void make_first_edge(site const l, site const r)
    {
        assert(endpoints_.empty());
        pedge const e = add_edge(l, r, nv);
        pendpoint const ll = insert_endpoint(nep, l, r, e);
        if (less_(l->x, r->x))  {
            pendpoint const rr = insert_endpoint(nep, r, l, e);
            assert(std::next(ll) == rr);
        }
    }

    void begin_cell(site const s, point const & _site)
    {
        pendpoint l = endpoints_.lower_bound(_site);
        pendpoint r = l;
        while (r != nep) {
            if (less_(_site, r->first)) {
                break;
            }
            assert(l->second.ev == r->second.ev); // if fires, then there is problem with precision
            ++r;
        }
        if (l == r) {
            if (l == nep) { // append to the rightmost endpoint
                --l;
                site const c = l->first.r;
                pedge const e = add_edge(c, s, nv);
                r = insert_endpoint(nep, c, s, e);
                if (less_(c->x, _site.x))  {
                    pendpoint const rr = insert_endpoint(nep, s, c, e);
                    assert(std::next(r) == rr);
                }
            } else if (l == std::begin(endpoints_)) { // prepend to the leftmost endpoint
                site const c = r->first.l;
                pedge const e = add_edge(s, c, nv);
                pendpoint const ll = insert_endpoint(r, c, s, e);
                l = insert_endpoint(r, s, c, e);
                assert(std::next(ll) == l);
            } else { // insert in the middle of the beachline (hottest branch in general case)
                --l;
                site const c = l->first.r;
                assert(c == r->first.l);
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
            if (endpoint_.second != nev) {
                assert(less_(_site.x, endpoint_.second.ev->first.x()));
                disable_event(endpoint_.second);
            }
            auto vertex_ = make_vertex(_site, *endpoint_.first.l, *endpoint_.first.r);
            assert(!!vertex_);
            assert(events_.find(*vertex_) == nev);
            pvertex const v = vertices_.insert(nv, std::move(*vertex_));
            trunc_edge(*endpoint_.first.e, v);
            pedge const le = add_edge(endpoint_.first.l, s, v);
            pedge const re = add_edge(s, endpoint_.first.r, v);
            pendpoint const ep = insert_endpoint(r, s, endpoint_.first.r, re);
            assert(std::next(ep) == r);
            endpoints_.erase(std::exchange(l, insert_endpoint(ep, endpoint_.first.l, s, le)));
            assert(std::next(l) == ep);
            if (l != std::begin(endpoints_)) {
                check_event(std::prev(l), l);
            }
            if (r != nep) {
                check_event(ep, r);
            }
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
    value_type
    angle(endpoint const & ep)
    {
        return angle(*ep.l, *ep.r);
    }

    std::pair< pendpoint, pendpoint >
    endpoint_range(pray const l, pray const r) const
    {
        assert(r != nray);
        assert(0 < std::distance(l, r));
        if (std::next(l) == r) {
            assert(std::next(*l) == *r);
            return {*l, *r};
        } else {
            auto const angle_less = [] (pendpoint const ll, pendpoint const rr) -> bool
            {
                return angle(ll->first) < angle(rr->first);
            };
            auto const lr = std::minmax_element(l, std::next(r), angle_less);
            return {*lr.first, *lr.second};
        }
    }

    bool check_endpoint_range(pevent const ev, pendpoint l, pendpoint const r) const
    {
        if (r == nep) {
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

    void finish_cells(pevent const ev,
                      vertex const & _vertex, bundle const & b,
                      site const l, site const r)
    {
        auto lr = endpoint_range(b.first, b.second);
        // All the edges from [*lr.first; *r.second]->first.e can be stored near the associate vertex if needed
        assert(check_endpoint_range(ev, lr.first, lr.second));
        pvertex const v = vertices_.insert(nv, _vertex);
        remove_event(ev, b);
        site const ll = lr.first->first.l;
        site const rr = lr.second->first.r;
        ++lr.second;
        do {
            trunc_edge(*lr.first->first.e, v);
            endpoints_.erase(lr.first++);
        } while (lr.first != lr.second);
        if (l == r) {
            lr.first = insert_endpoint(lr.second, ll, rr, add_edge(ll, rr, v));
            if (lr.first != std::begin(endpoints_)) {
                check_event(std::prev(lr.first), lr.first);
            }
            if (lr.second != nep) {
                check_event(lr.first, lr.second);
            }
        } else {
            pendpoint const ep = insert_endpoint(lr.second, l, rr, add_edge(l, rr, v));
            lr.first = insert_endpoint(ep, ll, l, add_edge(ll, l, v));
            assert(std::next(lr.first) == ep);
            if (lr.first != std::begin(endpoints_)) {
                check_event(std::prev(lr.first), lr.first);
            }
            if (lr.second != nep) {
                check_event(ep, lr.second);
            }
        }
    }

    bool check_last_endpoints() const
    {
        for (auto const & ep : endpoints_) {
            edge const & e = *ep.first.e;
            if ((e.b != nv) && (e.e != nv)) {
                return false;
            }
            if (ep.second != nev) {
                return false;
            }
        }
        return true;
    }

public :

    template< typename iterator >
    void operator () (iterator l, iterator const r)
    {
        assert(std::is_sorted(l, r, less_));
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
            point const & site_ = *l;
            bool next = false;
            while (!events_.empty()) {
                pevent const ev = std::begin(events_);
                auto const & event_ = *ev;
                {
                    value_type const & x = event_.first.x();
                    if (less_(site_.x, x)) {
                        break;
                    } else if (!less_(x, site_.x)) {
                        value_type const & y = event_.first.y();
                        if (less_(site_.y, y)) {
                            break;
                        } else if (!less_(y, site_.y)) {
                            finish_cells(ev, event_.first, event_.second, l, r);
                            next = true;
                            break;
                        }
                    }
                }
                finish_cells(ev, event_.first, event_.second, l, l);
            }
            if (!next) {
                begin_cell(l, site_);
            }
        }
        while (!events_.empty()) {
            pevent const ev = std::begin(events_);
            auto const & event_ = *ev;
            finish_cells(ev, event_.first, event_.second, r, r);
        }
        assert(std::is_sorted(std::begin(vertices_), nv, less_)); // bigger then usual roundoff errors may lead to fire, though not mutters much
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
