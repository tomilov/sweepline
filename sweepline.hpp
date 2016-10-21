#pragma once

#include <algorithm>
#include <deque>
#include <functional>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#if 1
#include <iostream>
#include <istream>
#include <ostream>
#endif

#include <cassert>
#include <cmath>

// TODO: map< k, v > instead of set< struct { k; mutable v; } >
template< typename point_iterator, typename value_type >
struct sweepline
{

    using size_type = std::size_t;

    value_type const & eps;

    sweepline(value_type const && _eps) = delete; // lifetime of sweepline instance must be exceeded by a lifetime of eps

    explicit
    sweepline(value_type const & _eps)
        : eps(_eps)
    { ; }

    using site = point_iterator;

    using sites = std::vector< site >;
    using site_iterator = typename sites::iterator;

    // TODO: add references to edges and/or sites
    struct vertex
    {

        value_type x, y;
        mutable value_type R;

    };

    struct vertex_less
    {

        value_type const & eps_;

        bool
        operator () (vertex const & _lhs, vertex const & _rhs) const // lexicographically compare
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
        }

    };

    using vertices = std::set< vertex, vertex_less >;
    using vertex_iterator = typename vertices::iterator;

    struct edge // segment, ray or line
    {

        site_iterator l, r; // for unfinished edges l, r, e is clockwise oriented triple
        vertex_iterator b, e; // end iterator may be invalidated even at move operation

    };

    using edges = std::list< edge >;
    using edge_iterator = typename edges::iterator;

    // input
    sites sites_;
    // output
    vertex_less vertex_less_{eps};
    vertices vertices_{vertex_less_};
    edges edges_;
    // TODO: clone() (not too stratiforward how to implement
    // due to (cross-referenced) past-the-end iterators invalidation when whole container movied)

    vertex_iterator const vend = std::end(vertices_);
    edge_iterator const inf = std::end(edges_);

    static
    value_type // ordinate
    intersect(edge const & _edge,
              bool const _ccw,
              value_type _directrix,
              value_type const & _eps)
    {
        auto const & u = **(_ccw ? _edge.l : _edge.r);
        auto const & v = **(_ccw ? _edge.r : _edge.l);
        {
            bool const v_degenerated_ = !(v.x + _eps < _directrix);
            if (!(u.x + _eps < _directrix)) {
                if (v_degenerated_) {
                    assert(u.y + _eps < v.y); // u != v
                    return (u.y + v.y) / value_type(2);
                } else {
                    return u.y;
                }
            } else if (v_degenerated_) {
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
        auto const calc_c = [&] (auto const & w, value_type const & wd)
        {
            return (w.x * w.x + w.y * w.y - _directrix) / wd;
        };
        value_type uc = calc_c(u, ud);
        value_type vc = calc_c(v, vd);
        value_type b = vb - ub; // -b
        value_type c = vc - uc;
        if ((u.x + _eps < v.x) || (v.x + _eps < u.x)) {
            value_type a = (ud - vd) / (ud * vd);
            a += a;
            value_type D = b * b - (a + a) * c;
            assert(!(D < value_type(0)));
            using std::sqrt;
            return (b + sqrt(std::move(D))) / a;
        } else { // a ~= 0
            return c / b; // -c / b
        }
    } // y(x)

    static
    value_type // absciss
    intersect(site const & _focus,
              value_type const & y,
              value_type const & _directrix)
    {
        auto const & focus_ = *_focus;
        assert(focus_.x < _directrix);
        value_type d = focus_.x - _directrix;
        return (y * (y - (focus_.y + focus_.y)) + (focus_.x * focus_.x + focus_.y * focus_.y - _directrix * _directrix)) / (d + d); // x(y)
    }

private :

    struct site_less
    {

        bool
        operator () (site const & _lhs, site const & _rhs) const // lexicographically compare w/o tolerance
        {
            auto const & lhs_ = *_lhs;
            auto const & rhs_ = *_rhs;
            return std::tie(lhs_.x, lhs_.y) < std::tie(rhs_.x, rhs_.y);
        }

    };

    edge_iterator
    start_edge(site_iterator const l, site_iterator const r, vertex_iterator const v) // begin from vertex
    {
        return edges_.insert(std::cend(edges_), {l, r, v, vend});
    }

    edge_iterator
    make_edge(site_iterator l, site_iterator r)
    {
        auto const & L = **l;
        auto const & R = **r;
        if ((R.y < L.y) || (!(L.y < R.y) && !(L.x < R.x))) {
            std::swap(l, r);
        }
        return start_edge(l, r, vend);
    }

    void
    finish_edge(edge & _edge, vertex_iterator const v) const
    {
        assert(v != vend);
        assert(_edge.e != v);
        if (_edge.b == vend) {
            _edge.b = v;
            auto const & l = **_edge.l;
            auto const & r = **_edge.r;
            vertex const & p = *v;
            bool swap = false;
            if (r.x < l.x) {
                if (l.y < p.y) {
                    swap = true;
                }
            } else if (l.x < r.x) {
                if (p.y < r.y) {
                    swap = true;
                }
            } else {
                swap = true;
            }
            if (swap) {
                std::swap(_edge.l, _edge.r);
            }
        } else {
            assert(_edge.b != v);
            assert(_edge.e == vend);
            _edge.e = v;
        }
    }

    // forward declaration for cross-referencing
    struct arc;
    struct arc_less;

    using beach_line = std::set< arc, arc_less >;
    using arc_iterator = typename beach_line::iterator;

    struct event
    {

        // non-copyable

        vertex_iterator circumcenter_; // vertex to be created
        value_type x;
        mutable arc_iterator l, r;     // inclusive range of collapsing arcs to be removed

        event(vertex_iterator _circumcenter,
              value_type && _x,
              arc_iterator const & a)
            : circumcenter_(std::move(_circumcenter))
            , x(std::move(_x))
            , l(a)
            , r(a)
        { ; }

        value_type const & y() const { return circumcenter_->y; }

    };

    struct event_less
    {

        using is_transparent = void;

        value_type const & eps_;

        bool
        operator () (event const & _lhs, event const & _rhs) const
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y() + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y());
        }

        bool
        operator () (event const & _lhs, vertex const & _rhs) const
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y() + eps_;
            value_type const & xx = _rhs.x + _rhs.R;
            return std::tie(x, y) < std::tie(xx, _rhs.y);
        }

        bool
        operator () (vertex const & _lhs, event const & _rhs) const
        {
            value_type const & x = _lhs.x + _lhs.R + eps_;
            value_type const & y = _lhs.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y());
        }

    };

    using events = std::set< event, event_less >;
    using event_iterator = typename events::iterator;

    struct arc
    {

        site_iterator focus_;
        mutable edge_iterator l, r;
        mutable event_iterator event_;

        friend
        std::ostream &
        operator << (std::ostream & _out, arc const & _arc)
        {
            auto const & p = *_arc.focus_->p;
            return _out << p.x << ' ' << p.y;
        }

    };

    struct arc_less
    {

        using is_transparent = void;

        value_type const & eps_;
        edge_iterator const & inf_;
        value_type const & directrix_;
        //bool instant_circle_event_ = false;

        bool
        operator () (arc const & _lhs, arc const & _rhs) const
        {
            if ((_lhs.r == inf_) || (_rhs.l == inf_)) {
                return false;
            }
            if (_lhs.r == _rhs.l) {
                return true;
            }
            if ((_rhs.r == inf_) || (_lhs.l == inf_)) {
                return true;
            }
            if (_rhs.r == _lhs.l) {
                return false;
            }
            edge const & lhs_ = *_lhs.r;
            edge const & rhs_ = *_rhs.l;
            value_type const l = intersect(lhs_, (_lhs.focus_ == lhs_.l), directrix_, eps_);
            value_type const r = intersect(rhs_, (_rhs.focus_ == rhs_.r), directrix_, eps_);
            return l + eps_ < r;
        }

        bool
        operator () (arc const & _lhs, site const & _rhs) const
        {
            if (_lhs.r == inf_) {
                return false;
            }
            edge const & lhs_ = *_lhs.r;
            auto const & focus_ = **_lhs.focus_;
            auto const & point_ = *_rhs;
            if (focus_.x + eps_ < point_.x) {
                value_type const intersection_ = intersect(lhs_, (_lhs.focus_ == lhs_.l), point_.x, eps_);
                if (intersection_ + eps_ < point_.y) {
                    return true;
                } else if (point_.y + eps_ < intersection_) {
                    return false;
                } else {
                    //instant_circle_event_ = true;
                    return false;
                }
            } else {
                auto const & right_ = **((_lhs.focus_ == lhs_.l) ? lhs_.r : lhs_.l);
                if (right_.x + eps_ < focus_.x) {
                    return focus_.y + eps_ < point_.y;
                } else {
                    assert(!(focus_.x + eps_ < right_.x)); // equivalent
                    return right_.y + eps_ < point_.y; // degenerated right neighbour // always true for lexicographically ordered and properly spaced sites
                }
            }
        }

        bool
        operator () (site const & _lhs, arc const & _rhs) const
        {
            if (_rhs.l == inf_) {
                return false;
            }
            edge const & rhs_ = *_rhs.l;
            auto const & focus_ = **_rhs.focus_;
            auto const & point_ = *_lhs;
            if (focus_.x + eps_ < point_.x) {
                value_type const intersection_ = intersect(rhs_, (_rhs.focus_ == rhs_.r), point_.x, eps_);
                if (point_.y + eps_ < intersection_) {
                    return true;
                } else if (intersection_ + eps_ < point_.y) {
                    return false;
                } else {
                    //instant_circle_event_ = true;
                    return false;
                }
            } else {
                auto const & left_ = **((_rhs.focus_ == rhs_.r) ? rhs_.l : rhs_.r);
                if (left_.x + eps_ < focus_.x) {
                    return point_.y + eps_ < focus_.y;
                } else {
                    assert(!(focus_.x + eps_ < left_.x)); // equivalent
                    return point_.y + eps_ < left_.y; // degenerated left neighbour // always false for lexicographically ordered and properly spaced sites
                }
            }
        }

    };

    value_type directrix_;
    arc_less arc_less_{eps, inf, directrix_};
    beach_line beach_line_{arc_less_};
    typename beach_line::const_iterator blend = std::cend(beach_line_);
    event_less event_less_{eps};
    events events_{event_less_};
    event_iterator const noe = std::end(events_);

    friend
    std::ostream &
    operator << (std::ostream & _out, beach_line const & b)
    {
        for (arc const & arc_ : b) {
            _out << arc_ << std::endl;
        }
        return _out << std::endl;
    }

    void
    delete_event(event_iterator const e)
    {
        if (e == noe) {
            return;
        }
        event const & event_ = *e;
        {
            auto a = event_.l;
            assert(a->event_ == e);
            while (a != event_.r) {
                if (a->event_ == e) {
                    a->event_ = noe;
                }
                ++a;
            }
            assert(a->event_ == e);
            a->event_ = noe;
        }
        vertices_.erase(event_.circumcenter_);
        events_.erase(e);
    }

    void
    disable_event(arc_iterator const & a)
    {
        event_iterator & e = a->event_;
        if (e == noe) {
            return;
        }
        event const & event_ = *e;
        if (event_.l == event_.r) {
            assert(event_.l == a);
            vertices_.erase(event_.circumcenter_);
            events_.erase(e);
        } else if (event_.l == a) {
            ++event_.l;
        } else if (event_.r == a) {
            --event_.r;
        } else {
            // TODO: implement
            assert(false);
        }
        e = noe;
    }

    void
    check_event(arc_iterator const & l,
                arc_iterator const & a,
                arc_iterator const & r)
    {
        assert(a->event_ == noe);
        assert(&l->focus_ != &a->focus_);
        assert(&r->focus_ != &a->focus_);
        assert(&l->focus_ != &r->focus_);
        auto const & u = **l->focus_;
        auto const & v = **a->focus_;
        auto const & w = **r->focus_;
        value_type A = v.x - u.x;
        value_type B = v.y - u.y;
        value_type C = w.x - u.x;
        value_type D = w.y - u.y;
        value_type G = B * (w.x - v.x) - A * (w.y - v.y);
        if (!(eps * eps < G)) {
            // 1.) non-concave triple of points => circumcircle don't cross the sweep line
            // 2.) G is small: collinear points => edges never cross
            return;
        }
        G += G;
        value_type E = A * (u.x + v.x) + B * (u.y + v.y);
        value_type F = C * (u.x + w.x) + D * (u.y + w.y);
        value_type x = (B * F - D * E) / G;
        value_type y = (C * E - A * F) / G;
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
        assert(eps * eps * eps < R); // are points too close to each other?
        R *= (e + f + g);
        R = (e * f * g) / sqrt(std::move(R));
        // R - radius
        auto const pvertex = vertices_.insert({x, y, R});
        if (pvertex.second) {
            auto const ee = events_.emplace(pvertex.first, x + R, a);
            assert(ee.second);
            a->event_ = ee.first;
        } else {
            vertex const & vertex_ = *pvertex.first;
            auto pevent = events_.find(vertex_);
            assert(pevent != noe);
            if (R + eps < vertex_.R) {
                vertex_.R = R;
                delete_event(pevent);
                auto const ee = events_.emplace(pvertex.first, x + R, a);
                assert(ee.second);
                a->event_ = ee.first;
            } else if (vertex_.R + eps < R) {
                // nop
            } else { // equiv
                event const & event_ = *pevent;
                if (r == event_.l) {
                    event_.l = a;
                } else if (l == event_.r) {
                    event_.r = a;
                } else {
                    if (arc_less_(*r, *event_.l)) {
                        event_.l = a;
                    } else if (arc_less_(*event_.r, *l)) {
                        event_.r = a;
                    } else {
                        assert(false); // leftmost arcs created first (if any), then upper and lower ones
                    }
                }
                a->event_ = pevent;
            }
        }
    }

    void
    add_arc(site_iterator const & s)
    {
        assert(!beach_line_.empty());
        site const & site_ = *s;
        auto const & point_ = *site_;
        directrix_ = point_.x;
        auto const range = beach_line_.equal_range(site_);
        assert(!arc_less_(*range.first, site_));
        assert(!arc_less_(site_, *range.first));
        assert(range.first != range.second); // because beach line is not empty
        auto const second = std::next(range.first);
        arc const & arc_ = *range.first;
        edge_iterator const ll = arc_.l;
        edge_iterator const rr = arc_.r;
        site_iterator const f = arc_.focus_;
        delete_event(arc_.event_);
        beach_line_.erase(range.first);
        if (second == range.second) { // 1 arc
            auto const & focus_ = **f;
            if (focus_.x + eps < point_.x) {
                edge_iterator const e = make_edge(f, s);
                arc_iterator const l = beach_line_.insert(second, {f, ll, e,  noe});
                arc_iterator const a = beach_line_.insert(second, {s, e,  e,  noe});
                arc_iterator const r = beach_line_.insert(second, {f, e,  rr, noe});
                assert(arc_less_(*l, *a));
                assert(arc_less_(*a, *r));
                if (l != std::cbegin(beach_line_)) {
                    check_event(std::prev(l), l, a);
                }
                if (second != blend) {
                    check_event(a, r, second);
                }
            } else if (second == blend) { // degenerated arc (horizontal line)
                assert(!(point_.x + eps < focus_.x));
                // horizontal line
                assert(focus_.y + eps < point_.y);
                edge_iterator const e = make_edge(f, s);
                beach_line_.insert(second, {f, ll, e, noe});
                beach_line_.insert(second, {s, e, inf, noe});
            } else {
                assert(false);
            }
        } else if (std::next(second) == range.second) { // 2 arcs
            //arc const ra = *second;
            beach_line_.erase(second);
            // TODO: hit the edge between two arcs, or even hit the vertex
            assert(false);
        } else {
            assert(false);
        }
    }

    void
    remove_arcs(event const & _event)
    {
        assert(!beach_line_.empty());
        assert(_event.r != blend);
        arc_iterator const ll = std::prev(_event.l);
        arc_iterator const rr = std::next(_event.r);
        assert(rr != blend);
        for (auto a = ll; a != rr; ++a) { // finish all the edges which are between interstitial arcs
            assert(a->r != inf);
            assert((a == ll) || (&*a->event_ == &_event));
            finish_edge(*a->r, _event.circumcenter_);
        }
        beach_line_.erase(_event.l, rr); // then remove supported arcs
        assert(std::next(ll) == rr);
        ll->r = rr->l = start_edge(ll->focus_, rr->focus_, _event.circumcenter_);
        //if (ll->event_ != rr->event_) {
        disable_event(ll);
        disable_event(rr);
        if (ll != std::cbegin(beach_line_)) {
            check_event(std::prev(ll), ll, rr);
        }
        if (std::next(rr) != blend) {
            check_event(ll, rr, std::next(rr));
        }
        //}
    }

public :

    void
    operator () ()
    {
        if (sites_.empty()) {
            return;
        }
        auto const send = std::end(sites_);
        std::sort(std::begin(sites_), send, site_less{});
        auto s = std::begin(sites_);
        beach_line_.insert({s, inf, inf, noe});
        while (++s != send) {
            if (!events_.empty()) {
                auto const & focus_ = **s;
                do {
                    auto const e = events_.begin();
                    event const & event_ = *e;
                    if (focus_.x + eps < event_.x) {
                        break;
                    }
                    if (!(event_.x + eps < focus_.x)) {
                        value_type const & y = event_.y();
                        if (!(y + eps < focus_.y)) {
                            if (!(focus_.y + eps < y)) {
                                break; // if event and focus are equivalent points, then site should be processed first
                            }
                        }
                    }
                    remove_arcs(event_);
                    events_.erase(e);
                } while (!events_.empty());
            }
            add_arc(s);
        }
        while (!events_.empty()) {
            auto const e = events_.begin();
            remove_arcs(*e);
            events_.erase(e);
        }
    }

};
