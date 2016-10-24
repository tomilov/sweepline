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
            return operator * () < *_rhs;
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

        // cached y:
        std::size_t n;
        value_type y;

    };

    struct endpoint_less
    {

        value_type const & eps_;
        value_type const & directrix_;

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

    using endpoints = std::set< endpoint, endpoint_less >;
    using pendpoint = typename endpoints::iterator;

    value_type directrix;
    endpoint_less endpoint_less_{eps, directrix};
    endpoints endpoints_{endpoint_less_};
    pendpoint const noe = std::end(endpoints_);

    struct event
    {

        pendpoint l, r;

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

    };

    using events = std::map< pvertex, event, event_less >;

    event_less event_less_{eps};
    events events_{event_less_};

    void
    begin_cell(site const &, cell const &)
    {

    }

    void
    finish_cell(pvertex const, event const &)
    {

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
