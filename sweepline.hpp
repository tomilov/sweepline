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

template< typename point_iterator, typename point_type, typename value_type >
struct sweepline
{

    static_assert(std::is_same_v< decltype(std::declval< point_type >().x), decltype(std::declval< point_type >().y) >, "point_type should have x and y data members");

    using size_type = std::size_t;

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    struct point_less
    {

        bool operator () (point_type const & _lhs, point_type const & _rhs) const // lexicographically compare w/o tolerance
        {
            return std::tie(_lhs.x, _lhs.y) < std::tie(_rhs.x, _rhs.y);
        }

    };

    struct site
    {

        point_iterator p;

        bool operator < (site const & _rhs) const
        {
            return point_less{}(*p, *_rhs.p);
        }

    };

    using vertices = std::map< point_type, value_type, point_less >; // maps vertices (circumcenter) to radii of circumscribed circles
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

    vertices vertices_;
    pvertex const nov = std::end(vertices_);
    edges edges_;
    cells cells_;

private :

    struct endpoint
    {

        pcell l, r;

    };

    struct endpoint_less
    {

        value_type const & eps_;

        bool operator () (endpoint const & _lhs, endpoint const & _rhs) const
        {
            if (_lhs.r == _rhs.l) {
                return true;
            }
            if (_lhs.l == _rhs.r) {
                return false;
            }
            return false;
        }

    };

    using endpoints = std::set< endpoint, endpoint_less >;

    endpoint_less endpoint_less_{eps};
    endpoints endpoints_{endpoint_less_};

public :

    void
    operator () (point_iterator beg, point_iterator const end)
    {
        if (beg == end) {
            return;
        }
        do {
            if (!cells_.insert({site{beg}, {}}).second) {
                assert(false); // bad input
            }
        } while (++beg != end);
    }

};
