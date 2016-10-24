#pragma once

#include <algorithm>
#include <deque>
#include <functional>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <set>
#include <map>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <experimental/optional>

#include <cassert>
#include <cmath>

static std::size_t numerator = 0;
static std::size_t denominator = 0;

template< typename point_iterator,
          typename point_type, // lexicagraphically comparable (x then y)
          typename value_type >
struct sweepline
{

    static_assert(std::is_same_v< decltype(std::declval< point_type >().x), decltype(std::declval< point_type >().y) >, "point_type should have x and y data members of the same type");

    using size_type = std::size_t;

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct site
    {

        site(point_iterator _p)
            : p(std::move(_p))
        { ; }

        point_iterator
        value() const
        {
            return p;
        }

        point_iterator
        operator -> () const
        {
            return p;
        }

        point_type &
        operator * () const
        {
            return *p;
        }

        bool operator < (site const & _rhs) const
        {
            return *p < *_rhs;
        }

    private :

        point_iterator p;

    };

    struct vertex_less
    {

        value_type const & eps_;

        bool operator () (point_type const & _lhs, point_type const & _rhs) const // lexicographically compare w/ tolerance
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
        }

    };

    using vertices = std::map< point_type, value_type, vertex_less >; // mapping of vertices (circumcenter) to radii of circumscribed circles
    using pvertex = typename vertices::iterator;

    struct edge
    {

        site l, r;
        pvertex b, e;

        bool operator < (edge const & _rhs) const
        {
            assert(l < r);
            assert(_rhs.l < _rhs.r);
            return std::tie(l, r) < std::tie(_rhs.l, _rhs.r);
        }

    };

    using edges = std::set< edge >; // or vector?
    using pedge = typename edges::iterator;

    using cell = std::deque< pedge >; // front() - left, back() - right
    using cells = std::map< site, cell >;
    using pcell = typename cells::iterator;

    vertex_less vertex_less_{eps};
    vertices vertices_{vertex_less_};
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_;

private :

    struct endpoint
    {

        pcell l, r;

    };

    struct event_less
    {

        value_type const & eps_;

        bool operator () (pvertex const & _lhs, pvertex const & _rhs) const
        {
            auto const & lhs_ = *_lhs;
            auto const & rhs_ = *_rhs;
            return std::make_pair((lhs_.first.x + lhs_.second) + eps_, lhs_.first.y + eps_) < std::make_pair(rhs_.first.x + rhs_.second, rhs_.first.y);
        }

        using is_transparent = void;

        bool operator () (pvertex const & _lhs, site const & _rhs) const
        {
            auto const & lhs_ = *_lhs;
            value_type const & x = (lhs_.first.x + lhs_.second) + eps_;
            value_type const & y = lhs_.first.y + eps_;
            auto const & rhs_ = *_rhs;
            return std::tie(x, y) < std::tie(rhs_.x, rhs_.y);
        }

        bool operator () (site const & _lhs, pvertex const & _rhs) const
        {
            auto const & lhs_ = *_lhs;
            auto const & rhs_ = *_rhs;
            return std::make_pair(lhs_.x + eps_, lhs_.y + eps_) < std::make_pair(rhs_.first.x + rhs_.second, rhs_.first.y);
        }

        value_type const & directrix_;
        std::size_t const & n;

        bool operator () (endpoint const & _lhs, endpoint const & _rhs) const
        {
            if (_lhs.r == _rhs.l) {
                return true;
            }
            if (_lhs.l == _rhs.r) {
                return false;
            }
            assert(false); // TODO: implement
        }

    };

    using endpoints = std::map< endpoint, pvertex, event_less >;
    using pendpoint = typename endpoints::iterator;

    struct endpoint_range
    {

        pendpoint b, e;

        pendpoint begin() const { return b; }
        pendpoint end()   const { return e; }

    };

    using events = std::map< pvertex, endpoint_range, event_less >;

    value_type directrix;
    std::size_t n;
    event_less event_less_{eps, directrix, n};
    endpoints endpoints_{event_less_};
    events events_{event_less_};

    void
    begin_cell(site const &, cell const &)
    {

    }

    std::experimental::optional< typename vertices::value_type >
    make_vertex(point_type const & u,
                point_type const & v,
                point_type const & w) const
    {
        value_type A = v.x - u.x;
        value_type B = v.y - u.y;
        value_type C = w.x - u.x;
        value_type D = w.y - u.y;
        value_type G = B * (w.x - v.x) - A * (w.y - v.y);
        if (!(eps * eps < G)) {
            // 1.) non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return {};
        }
        G += G;
        value_type E = A * (u.x + v.x) + B * (u.y + v.y);
        value_type F = C * (u.x + w.x) + D * (u.y + w.y);
        point_type point_{(B * F - D * E) / G, (C * E - A * F) / G};
        // x, y - circumcenter
        using std::sqrt;
        auto const norm = [&] (auto const & ll, auto const & rr) -> value_type
        {
            value_type dx = rr.x - ll.x;
            value_type dy = rr.y - ll.y;
            return sqrt(dx * dx + dy * dy);
        };
        value_type e = norm(u, v);
        value_type f = norm(v, w);
        value_type g = norm(w, u);
        value_type R = (e + f - g) * (e + g - f) * (f + g - e);
        assert(eps * eps * eps < R);
        R = (e * f * g) / sqrt(std::move(R) * (e + f + g));
        // R - radius
        return {{std::move(point_), std::move(R)}};
    }

    void
    replace_event(pvertex & v, pvertex V) const
    {
        if (v == nov) ; else if (event_less_(V, v)) {
            // disable event (whatever it mean, I tired)
        } else {
            return;
        }
        v = std::move(V);
    }

    void
    check_event(pendpoint const l, pendpoint const r)
    {
        assert(std::next(l) == r);
        assert(l->first.r == r->first.l);
        if (auto vertex_ = make_vertex(*l->first.l->first, *l->first.r->first, *r->first.r->first)) {
            auto const v = vertices_.insert(std::move(*vertex_));
            replace_event(l->second, v.first);
            replace_event(r->second, v.first);
            if (v.second) {
                if (!events_.insert({v.first, {l, std::next(r)}}).second) {
                    assert(false);
                }
            } else {
                // implement
                assert(false);
            }
        }
    }

    void
    finish_cell(pvertex const & _circumcenter, endpoint_range _range)
    {
        assert(1 < endpoints_.size()); // endpoints are removed at least by two
        assert(1 < std::distance(_range.b, _range.e));
        for (auto const & endpoint_ : _range) {
            // finish right edge of *endpoint_.first.l at _circumcenter vertex
            (void)_circumcenter;
            (void)endpoint_;
            assert(false); // TODO: implement
        }
        endpoint endpoint_{_range.b->first.l, std::prev(_range.e)->first.r};
        endpoints_.erase(_range.b, _range.e);
        // make new edge
        // add new edge to endpoint_.r->second.push_front() and to endpoint_.l->second.push_back()
        assert(endpoints_.find(endpoint_) == std::end(endpoints_));
        _range.b = endpoints_.insert(_range.e, {std::move(endpoint_), nov});
        // disable events for endpoint_
        if (_range.b != std::begin(endpoints_)) {
            check_event(std::prev(_range.b), _range.b);
        }
        if (_range.e != std::end(endpoints_)) {
            check_event(_range.b, _range.e);
        }
    }

public :

    bool
    operator () (point_iterator beg, point_iterator const end)
    {
        if (beg == end) {
            return true;
        }
        do {
            if (!cells_.insert({beg, {}}).second) {
                return false;
            }
        } while (++beg != end);
        for (auto & cell_ : cells_) {
            while (!events_.empty()) {
                auto const e = std::cbegin(events_);
                auto const & event_ = *e;
                if (!event_less_(event_.first, cell_.first)) {
                    break;
                }
                finish_cell(event_.first, event_.second);
                events_.erase(e);
            }
            begin_cell(cell_.first, cell_.second);
        }
        while (!events_.empty()) {
            auto const e = std::cbegin(events_);
            auto const & event_ = *e;
            finish_cell(event_.first, event_.second);
            events_.erase(e);
        }
        return true;
    }

};
