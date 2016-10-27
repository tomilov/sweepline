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

    static_assert(std::is_same< decltype(std::declval< point_type >().x), decltype(std::declval< point_type >().y) >::value, "point_type should have x and y data members of the same type");

    using size_type = std::size_t;

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct site_less
    {

        bool operator () (point_iterator const l, point_iterator const r) const
        {
            return *l < *r;
        }

    };

    struct vertex // denote circumscribed circle
    {

        point_type p; // circumcenter
        value_type R; // radius

        value_type x() const { return p.x + R; }
        value_type const & y() const { return p.y; }

    };

    struct vertex_less
    {

        value_type const & eps_;

        bool operator () (vertex const & l, vertex const & r) const // lexicographically compare w/ tolerance
        {
            value_type const & x = l.p.x + eps_;
            value_type const & y = l.p.y + eps_;
            return std::tie(x, y) < std::tie(r.p.x, r.p.y);
        }

    };

    using vertices = std::set< vertex, vertex_less >;
    using pvertex = typename vertices::iterator;

    struct edge // ((b, e), (l, r)) is ccw
    {

        point_iterator l, r;
        pvertex b, e;

    };

    using edges = std::list< edge >;
    using pedge = typename edges::iterator;

    using cells = std::map< point_iterator, std::deque< pedge >, site_less >;
    using pcell = typename cells::iterator;

    vertices vertices_{vertex_less{eps}};
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_;

private :

    struct endpoint
    {

        pcell l, r;
        pedge e;

    };

    value_type directrix = std::numeric_limits< value_type >::quiet_NaN();

    struct local_insert_hint
    {

        endpoint const * l = nullptr;
        endpoint const * m = nullptr;
        endpoint const * r = nullptr;

        void clear() { l = m = r = nullptr; }

    };

    struct endpoint_less
    {

        value_type const & eps_;

        value_type const & directrix_;

        value_type // ordinate
        intersect(point_type const & u,
                  point_type const & v,
                  value_type _directrix) const
        {
            {
                bool const degenerated_ = !(v.x + eps_ < _directrix);
                if (!(u.x + eps_ < _directrix)) {
                    if (degenerated_) {
                        assert(u.y + eps_ < v.y); // u != v
                        return (u.y + v.y) / value_type(2);
                    } else {
                        return u.y;
                    }
                } else if (degenerated_) {
                    return v.y;
                }
            }
            value_type ud = u.x - _directrix;
            value_type vd = v.x - _directrix;
            value_type ub = u.y / ud; // -b
            value_type vb = v.y / vd; // -b
            ud += ud;
            vd += vd;
            _directrix *= _directrix;
            auto const calc_c = [&] (point_type const & w, value_type const & wd)
            {
                return (w.x * w.x + w.y * w.y - _directrix) / wd;
            };
            value_type uc = calc_c(u, ud);
            value_type vc = calc_c(v, vd);
            value_type b = vb - ub; // -b
            value_type c = vc - uc;
            if ((u.x + eps_ < v.x) || (v.x + eps_ < u.x)) {
                value_type a = (ud - vd) / (ud * vd);
                a += a;
                value_type D = b * b - (a + a) * c;
                assert(!(D < value_type(0)));
                using std::sqrt;
                return (b + sqrt(D)) / a;
            } else { // a ~= 0
                return c / b; // -c / b
            }
        }

        local_insert_hint const & hint_;

        bool operator () (endpoint const & l, endpoint const & r) const
        {
            if (&l == &r) {
                return false;
            }
            if (hint_.m) {
                if (hint_.m == &l) {
                    if (hint_.l == &r) {
                        return false;
                    }
                    if (!hint_.r) {
                        return false;
                    }
                    if (hint_.r == &r) {
                        return true;
                    }
                } else if (hint_.m == &r) {
                    if (hint_.r == &l) {
                        return false;
                    }
                    if (!hint_.l) {
                        return false;
                    }
                    if (hint_.l == &r) {
                        return true;
                    }
                }
            }
            if (l.r == r.l) {
                return true;
            }
            if (l.l == r.r) {
                return false;
            }
            return intersect(*l.l->first, *l.r->first, directrix_) + eps_ < intersect(*r.l->first, *r.r->first, directrix_);
        }

        using is_transparent = void;

        bool operator () (vertex const & l, endpoint const & r) const
        { // need to replace with linear if possible
            return l.y() + eps_ < intersect(*r.l->first, *r.r->first, l.x());
        }

        bool operator () (endpoint const & l, vertex const & r) const
        { // need to replace with linear if possible
            return intersect(*l.l->first, *l.r->first, r.x()) + eps_ < r.y();
        }

    };

    using endpoints = std::map< endpoint, pvertex, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    local_insert_hint hint;
    endpoint_less const endpoint_less_{eps, directrix, hint};
    endpoints endpoints_{endpoint_less_};
    pendpoint const noep = std::end(endpoints_);

    struct event_less
    {

        value_type const & eps_;

        bool operator () (pvertex const & l, pvertex const & r) const
        { // uncomment if not all the events in equal range are for equally the same vertex
            vertex const & ll = *l;
            value_type const & x = ll.x() + eps_;
            value_type const & y = ll.y() + eps_;
            //value_type const & px = lhs_.p.x + eps_;
            vertex const & rr = *r;
            value_type const & vx = rr.x();
            return std::tie(x, y/*, px*/) < std::tie(vx, rr.y()/*, rhs_.p.x*/);
        }

    };

    using events = std::map< pvertex, pendpoint, event_less >;
    using pevent = typename events::const_iterator;

    event_less const event_less_{eps};
    events events_{event_less_};

    void
    begin_cell(point_iterator const p)
    {
        assert(cells_.find(p) == cells_.end());
        pcell const c = cells_.insert({p, {}}).first;
        (void)c; // TODO: implement
    }

    value_type
    circumradius(point_type const & u,
                 point_type const & v,
                 point_type const & w) const
    {
        using std::hypot;
        value_type e = hypot(v.x - u.x, v.y - u.y);
        value_type f = hypot(v.x - w.x, v.y - w.y);
        value_type g = hypot(w.x - u.x, w.y - u.y);
        value_type V = (e + f - g) * (e + g - f) * (f + g - e);
        assert(eps < V); // triangle inequality
        using std::sqrt;
        return (e * f * g) / sqrt(V * (e + f + g));
    }

    pvertex
    make_vertex(point_type const & u,
                point_type const & v,
                point_type const & w)
    {
        value_type A = v.x - u.x;
        value_type B = v.y - u.y;
        value_type C = w.x - u.x;
        value_type D = w.y - u.y;
        value_type G = B * (w.x - v.x) - A * (w.y - v.y);
        if (!(eps * eps < G)) {
            // 1.) non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return nov;
        }
        value_type E = A * (u.x + v.x) + B * (u.y + v.y);
        value_type F = C * (u.x + w.x) + D * (u.y + w.y);
        G += G;
        // circumcenter:
        value_type x = (B * F - D * E) / G;
        value_type y = (C * E - A * F) / G;
        return vertices_.insert({{x, y}, circumradius(u, v, w)}).first;
    }

    std::pair< pendpoint, pendpoint >
    endpoint_range(pendpoint const ep, pvertex const v)
    {
        if (ep == noep) {
            return endpoints_.equal_range(*v);
        } else {
            assert(ep != std::begin(endpoints_));
            assert(ep != noep);
            assert(std::prev(ep)->second == v);
            assert(ep->second == v);
            return {std::prev(ep), std::next(ep)};
        }
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

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        auto & ll = *l;
        auto & rr = *r;
        assert(ll.first.r == rr.first.l); // small redundancy (x2 memory to point cells from beach line) - necessary evil
        pvertex const v = make_vertex(*ll.first.l->first, *ll.first.r->first, *rr.first.r->first);
        if (v != nov) {
            if (ll.second != nov) {
                assert(rr.second == nov);
                if (event_less_(v, ll.second)) {
                    delete_event(ll.second);
                } else {
                    vertices_.erase(v);
                    return;
                }
            } else if (rr.second != nov) {
                assert(ll.second == nov);
                if (event_less_(v, rr.second)) {
                    delete_event(rr.second);
                } else {
                    vertices_.erase(v);
                    return;
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
    }

    void
    finish_edge(edge & _edge, pvertex const v) const
    {
        assert(v != nov);
        assert(_edge.e == nov);
        if (_edge.b == nov) {
            _edge.b = v;
            // adjust orinetation if needed:
            point_type const & l = *_edge.l;
            point_type const & r = *_edge.r;
            point_type const & p = v->p;
            if (r.x < l.x) {
                if (p.y < l.y) {
                    return;
                }
            } else if (l.x < r.x) {
                if (r.y < p.y) {
                    return;
                }
            } else {
                assert(!(r.y < l.y));
            }
            std::swap(_edge.l, _edge.r);
        } else {
            assert(_edge.b != v);
            _edge.e = v;
        }
    }

    void
    finish_edges(pevent const e)
    {
        auto const & event_ = *e;
        auto lr = endpoint_range(event_.second, event_.first);
        events_.erase(e);
        pcell const lcell = lr.first->first.l;
        pcell const rcell = std::prev(lr.second)->first.r;
        for (pendpoint ep = lr.first; ep != lr.second; ++ep) {
            finish_edge(*ep->first.e, event_.first);
        }
        endpoints_.erase(lr.first, lr.second);
        pedge const edge_ = edges_.insert(std::cend(edges_), {lcell->first, rcell->first, event_.first, nov});
        lr.first = endpoints_.insert(lr.second, {{lcell, rcell, edge_}, nov});
        lcell->second.push_front(edge_); // ccw
        rcell->second.push_back(edge_);
        if (lr.first != std::begin(endpoints_)) {
            check_event(std::prev(lr.first), lr.first);
        }
        if (lr.second != noep) {
            check_event(lr.first, lr.second);
        }
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
    operator () (point_iterator const beg, point_iterator const end)
    {
        for (auto p = beg; p != end; ++p) {
            while (!events_.empty()) {
                pevent const e = std::cbegin(events_);
                if (!prior(*e->first, *p)) {
                    break;
                }
                finish_edges(e);
            }
            begin_cell(p);
        }
        while (!events_.empty()) {
            finish_edges(std::cbegin(events_));
        }
    }

};
