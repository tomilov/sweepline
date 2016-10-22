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

    // TODO: add references to edges and/or sites
    struct point
    {

        value_type x, y;

    };

    struct point_less
    {

        value_type const & eps_;

        bool
        operator () (point const & _lhs, point const & _rhs) const // lexicographically compare w/ tolerance
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y);
        }

    };

    using vertices = std::map< point, value_type, point_less >;
    using vertex = typename vertices::value_type;
    using vertex_iterator = typename vertices::iterator;

    struct edge // segment, ray or line
    {

        site l, r; // for unfinished edges l, r, e is clockwise oriented triple
        vertex_iterator b, e; // end iterator may be invalidated even at move operation

    };

    using edges = std::list< edge >;
    using edge_iterator = typename edges::iterator;

    point_less point_less_{eps};
    vertices vertices_{point_less_};
    vertex_iterator const vend = std::end(vertices_);

    edges edges_;
    edge_iterator const inf = std::end(edges_);

    struct side
    {

        site const & l;
        site const & r;

    };

    static
    value_type // ordinate
    intersect(side const & _side,
              value_type _directrix,
              value_type const & _eps)
    {
        auto const & u = *_side.l;
        auto const & v = *_side.r;
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

    edge_iterator
    start_edge(site const l, site const r, vertex_iterator const v) // begin from vertex
    {
        return edges_.insert(std::cend(edges_), {l, r, v, vend});
    }

    edge_iterator
    make_edge(site const l, site const r)
    {
        return start_edge(l, r, vend);
    }

    void
    finish_edge(edge & _edge, vertex_iterator const v) const
    {
        assert(v != vend);
        assert(_edge.e != v);
        if (_edge.b == vend) {
            _edge.b = v;
            auto const & l = *_edge.l;
            auto const & r = *_edge.r;
            point const & p = v->first;
            if (r.x < l.x) {
                if (p.y < l.y) {
                    return;
                }
            } else if (l.x < r.x) {
                if (r.y < p.y) {
                    return;
                }
            } else {
                if (r.y < l.y) {
                    return; // unrecheable
                }
            }
            std::swap(_edge.l, _edge.r);
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

    struct event_location
    {

        // non-copyable

        vertex_iterator circumcenter_; // vertex to be created
        value_type x;

        event_location(vertex_iterator _circumcenter,
              value_type && _x)
            : circumcenter_(std::move(_circumcenter))
            , x(std::move(_x))
        { ; }

        value_type const & y() const { return circumcenter_->first.y; }

    };

    struct arc_range
    {

        arc_iterator l, r = l;     // inclusive range of collapsing arcs to be removed

    };

    struct event_less
    {

        using is_transparent = void;

        value_type const & eps_;

        bool
        operator () (event_location const & _lhs, event_location const & _rhs) const
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y() + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y());
        }

        bool
        operator () (event_location const & _lhs, vertex const & _rhs) const
        {
            value_type const & x = _lhs.x + eps_;
            value_type const & y = _lhs.y() + eps_;
            value_type const & xx = _rhs.first.x + _rhs.second;
            return std::tie(x, y) < std::tie(xx, _rhs.first.y);
        }

        bool
        operator () (vertex const & _lhs, event_location const & _rhs) const
        {
            value_type const & x = _lhs.first.x + _lhs.second + eps_;
            value_type const & y = _lhs.first.y + eps_;
            return std::tie(x, y) < std::tie(_rhs.x, _rhs.y());
        }

    };

    using events = std::map< event_location, arc_range, event_less >;
    using event = typename events::value_type;
    using event_iterator = typename events::iterator;

    struct arc
    {

        site focus_;
        mutable edge_iterator l, r;
        mutable event_iterator event_;

    };

    static
    side
    left(arc const & _arc)
    {
        edge const & edge_ = *_arc.r;
        if (_arc.focus_ == edge_.l) {
            return {edge_.l, edge_.r};
        } else {
            return {edge_.r, edge_.l};
        }
    }

    static
    side
    right(arc const & _arc)
    {
        edge const & edge_ = *_arc.l;
        if (_arc.focus_ == edge_.r) {
            return {edge_.l, edge_.r};
        } else {
            return {edge_.r, edge_.l};
        }
    }

    struct arc_less
    {

        using is_transparent = void;

        value_type const & eps_;
        value_type const & directrix_;
        edge_iterator const & inf_;
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
            return intersect(left(_lhs), directrix_, eps_) + eps_ < intersect(right(_rhs), directrix_, eps_);
        }

        bool
        operator () (arc const & _lhs, site const & _rhs) const
        {
            if (_lhs.r == inf_) {
                return false;
            }
            auto const & focus_ = *_lhs.focus_;
            auto const & point_ = *_rhs;
            side const left_ = left(_lhs);
            if (focus_.x + eps_ < point_.x) {
                value_type const intersection_ = intersect(left_, point_.x, eps_);
                if (intersection_ + eps_ < point_.y) {
                    return true;
                } else if (point_.y + eps_ < intersection_) {
                    return false;
                } else {
                    //instant_circle_event_ = true;
                    return false;
                }
            } else {
                auto const & right_ = *left_.r;
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
            auto const & focus_ = *_rhs.focus_;
            auto const & point_ = *_lhs;
            side const right_ = right(_rhs);
            if (focus_.x + eps_ < point_.x) {
                value_type const intersection_ = intersect(right_, point_.x, eps_);
                if (point_.y + eps_ < intersection_) {
                    return true;
                } else if (intersection_ + eps_ < point_.y) {
                    return false;
                } else {
                    //instant_circle_event_ = true;
                    return false;
                }
            } else {
                auto const & left_ = *right_.l;
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

    arc_less const arc_less_{eps, directrix_, inf};
    beach_line beach_line_{arc_less_};
    arc_iterator const brink = std::cend(beach_line_);

    event_less event_less_{eps};
    events events_{event_less_};
    event_iterator const noe = std::end(events_);

    void
    delete_event(event_iterator const e)
    {
        if (e == noe) {
            return;
        }
        arc_range const & arc_range_ = e->second;
        {
            auto a = arc_range_.l;
            assert(a->event_ == e);
            while (a != arc_range_.r) {
                if (a->event_ == e) {
                    a->event_ = noe;
                }
                ++a;
            }
            assert(a->event_ == e);
            a->event_ = noe;
        }
        vertices_.erase(e->first.circumcenter_);
        events_.erase(e);
    }

    void
    disable_event(arc_iterator const & a)
    {
        event_iterator & e = a->event_;
        if (e == noe) {
            return;
        }
        arc_range & arc_range_ = e->second;
        if (arc_range_.l == arc_range_.r) {
            assert(arc_range_.l == a);
            vertices_.erase(e->first.circumcenter_);
            events_.erase(e);
        } else if (arc_range_.l == a) {
            ++arc_range_.l;
        } else if (arc_range_.r == a) {
            --arc_range_.r;
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
        auto const & u = *l->focus_;
        auto const & v = *a->focus_;
        auto const & w = *r->focus_;
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
        auto const pvertex = vertices_.insert({{x, y}, R});
        if (pvertex.second) {
            auto const ee = events_.insert({{pvertex.first, x + R}, {a}});
            assert(ee.second);
            a->event_ = ee.first;
        } else {
            vertex & vertex_ = *pvertex.first;
            auto pevent = events_.find(vertex_);
            assert(pevent != noe);
            if (R + eps < vertex_.second) {
                vertex_.second = R;
                delete_event(pevent);
                auto const ee = events_.insert({{pvertex.first, x + R}, {a}});
                assert(ee.second);
                a->event_ = ee.first;
            } else if (vertex_.second + eps < R) {
                // nop
            } else { // equiv
                arc_range & arc_range_ = pevent->second;
                if (r == arc_range_.l) {
                    arc_range_.l = a;
                } else if (l == arc_range_.r) {
                    arc_range_.r = a;
                } else {
                    if (arc_less_(*r, *arc_range_.l)) {
                        arc_range_.l = a;
                    } else if (arc_less_(*arc_range_.r, *l)) {
                        arc_range_.r = a;
                    } else {
                        assert(false); // leftmost arcs created first (if any), then upper and lower ones
                    }
                }
                a->event_ = pevent;
            }
        }
    }

    void
    add_arc(site const & _site)
    {
        assert(!beach_line_.empty());
        auto const & point_ = *_site;
        directrix_ = point_.x;
        auto const range = beach_line_.equal_range(_site);
        assert(!arc_less_(*range.first, _site));
        assert(!arc_less_(_site, *range.first));
        assert(range.first != range.second); // because beach line is not empty
        auto const second = std::next(range.first);
        arc const & arc_ = *range.first;
        edge_iterator const ll = arc_.l;
        edge_iterator const rr = arc_.r;
        site const f = arc_.focus_;
        delete_event(arc_.event_);
        beach_line_.erase(range.first);
        if (second == range.second) { // 1 arc
            auto const & focus_ = *f;
            if (focus_.x + eps < point_.x) {
                edge_iterator const e = make_edge(f, _site);
                arc_iterator const l = beach_line_.insert(second, {f, ll, e,  noe});
                arc_iterator const a = beach_line_.insert(second, {_site, e,  e,  noe});
                arc_iterator const r = beach_line_.insert(second, {f, e,  rr, noe});
                assert(arc_less_(*l, *a));
                assert(arc_less_(*a, *r));
                if (l != std::cbegin(beach_line_)) {
                    check_event(std::prev(l), l, a);
                }
                if (second != brink) {
                    check_event(a, r, second);
                }
            } else if (second == brink) { // degenerated arc (horizontal line)
                assert(!(point_.x + eps < focus_.x));
                // horizontal line
                assert(focus_.y + eps < point_.y);
                edge_iterator const e = make_edge(f, _site);
                beach_line_.insert(second, {f, ll, e, noe});
                beach_line_.insert(second, {_site, e, inf, noe});
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
    remove_arcs(event & _event)
    {
        assert(!beach_line_.empty());
        assert(_event.second.r != brink);
        arc_iterator const ll = std::prev(_event.second.l);
        arc_iterator const rr = std::next(_event.second.r);
        assert(rr != brink);
        for (auto a = ll; a != rr; ++a) { // finish all the edges which are between interstitial arcs
            assert(a->r != inf);
            assert((a == ll) || (&*a->event_ == &_event));
            finish_edge(*a->r, _event.first.circumcenter_);
        }
        beach_line_.erase(_event.second.l, rr); // then remove supported arcs
        assert(std::next(ll) == rr);
        ll->r = rr->l = start_edge(ll->focus_, rr->focus_, _event.first.circumcenter_);
        disable_event(ll);
        disable_event(rr);
        if (ll != std::cbegin(beach_line_)) {
            check_event(std::prev(ll), ll, rr);
        }
        if (std::next(rr) != brink) {
            check_event(ll, rr, std::next(rr));
        }
    }

public :

    void
    operator () (point_iterator beg, point_iterator const end)
    {
        if (beg == end) {
            return;
        }
        std::deque< site > sites_;
        do {
            sites_.push_back(beg);
        } while (++beg != end);
        auto const site_less_ = [&] (site const & _lhs, site const & _rhs) -> bool
        {
            auto const & lhs_ = *_lhs;
            auto const & rhs_ = *_rhs;
            return std::tie(lhs_.x, lhs_.y) < std::tie(rhs_.x, rhs_.y); // lexicographically compare w/o tolerance
        };
        std::sort(std::begin(sites_), std::end(sites_), site_less_);
        beach_line_.insert({sites_.front(), inf, inf, noe});
        sites_.pop_front();
        for (site const & site_ : sites_) {
            if (!events_.empty()) {
                auto const & focus_ = *site_;
                do {
                    auto const e = events_.begin();
                    event & event_ = *e;
                    if (focus_.x + eps < event_.first.x) {
                        break;
                    }
                    if (!(event_.first.x + eps < focus_.x)) {
                        value_type const & y = event_.first.y();
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
            add_arc(site_);
        }
        while (!events_.empty()) {
            auto const e = events_.begin();
            remove_arcs(*e);
            events_.erase(e);
        }
    }

};
