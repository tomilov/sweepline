#pragma once

#include <type_traits>
#include <utility>
#include <iterator>
#include <memory>

#include <cassert>

namespace rb_tree
{

template< typename K, typename V >
struct pair { K k; V v; };

enum class color : bool { red = false, black = true };

struct node_base;

using base_pointer = node_base *;

struct node_base
{

    base_pointer p;
    base_pointer l;
    base_pointer r;
    color c;

    node_base() = default;

    node_base(base_pointer self) noexcept
        : p{nullptr}
        , l{self}
        , r{self}
        , c{color::red}
    { ; }

};

inline base_pointer minimum(base_pointer x) noexcept { while (x->l) { x = x->l; } return x; }
inline base_pointer maximum(base_pointer x) noexcept { while (x->r) { x = x->r; } return x; }

inline
base_pointer
increment(base_pointer x) noexcept
{
    if (x->r) {
        x = minimum(x->r);
    } else {
        base_pointer y = x->p;
        while (x == y->r) {
            x = y;
            y = y->p;
        }
        if (x->r != y) {
            x = y;
        }
    }
    return x;
}

inline
base_pointer
decrement(base_pointer x) noexcept
{
    if ((x->c == color::red) && (x->p->p == x)) {
        x = x->r;
    } else if (x->l) {
        x = maximum(x->l);
    } else {
        base_pointer y = x->p;
        while (x == y->l) {
            x = y;
            y = y->p;
        }
        x = y;
    }
    return x;
}

inline
void rotate_left(const base_pointer x, base_pointer & root) noexcept
{
    const base_pointer y = x->r;
    x->r = y->l;
    if (y->l) {
        y->l->p = x;
    }
    y->p = x->p;
    if (x == root) {
        root = y;
    } else if (x == x->p->l) {
        x->p->l = y;
    } else {
        x->p->r = y;
    }
    y->l = x;
    x->p = y;
}

inline
void rotate_right(const base_pointer x, base_pointer & root) noexcept
{
    const base_pointer y = x->l;
    x->l = y->r;
    if (y->r) {
        y->r->p = x;
    }
    y->p = x->p;
    if (x == root) {
        root = y;
    } else if (x == x->p->r) {
        x->p->r = y;
    } else {
        x->p->l = y;
    }
    y->r = x;
    x->p = y;
}

inline
void insert_and_rebalance(const bool insert_left,
                          base_pointer x, const base_pointer p,
                          node_base & h) noexcept
{
    base_pointer & root = h.p;
    x->p = p;
    x->l = nullptr;
    x->r = nullptr;
    x->c = color::red;
    if (insert_left) {
        p->l = x;
        if (p == &h) {
            h.p = x;
            h.r = x;
        } else if (p == h.l) {
            h.l = x;
        }
    } else {
        p->r = x;
        if (p == h.r) {
            h.r = x;
        }
    }
    while ((x != root) && (x->p->c == color::red)) {
        const base_pointer xpp = x->p->p;
        if (x->p == xpp->l) {
            const base_pointer y = xpp->r;
            if (y && (y->c == color::red)) {
                x->p->c = color::black;
                y->c = color::black;
                xpp->c = color::red;
                x = xpp;
            } else {
                if (x == x->p->r) {
                    x = x->p;
                    rotate_left(x, root);
                }
                x->p->c = color::black;
                xpp->c = color::red;
                rotate_right(xpp, root);
            }
        } else {
            const base_pointer y = xpp->l;
            if (y && y->c == color::red) {
                x->p->c = color::black;
                y->c = color::black;
                xpp->c = color::red;
                x = xpp;
            } else {
                if (x == x->p->l) {
                    x = x->p;
                    rotate_right(x, root);
                }
                x->p->c = color::black;
                xpp->c = color::red;
                rotate_left(xpp, root);
            }
        }
    }
    root->c = color::black;
}

inline
base_pointer
rebalance_for_erase(const base_pointer z, node_base & h) noexcept
{
    base_pointer & root = h.p;
    base_pointer & leftmost = h.l;
    base_pointer & rightmost = h.r;
    base_pointer y = z;
    base_pointer x = nullptr;
    base_pointer xp = nullptr;
    if (y->l) {
        if (y->r) {
            y = minimum(y->r);
            x = y->r;
        } else {
            x = y->l;
        }
    } else {
        x = y->r;
    }
    if (y != z) {
        z->l->p = y;
        y->l = z->l;
        if (y != z->r) {
            xp = y->p;
            if (x) {
                x->p = y->p;
            }
            y->p->l = x;
            y->r = z->r;
            z->r->p = y;
        } else {
            xp = y;
        }
        if (root == z) {
            root = y;
        } else if (z->p->l == z) {
            z->p->l = y;
        } else {
            z->p->r = y;
        }
        y->p = z->p;
        std::swap(y->c, z->c);
        y = z;
    } else {
        xp = y->p;
        if (x) {
            x->p = y->p;
        }
        if (root == z) {
            root = x;
        } else {
            if (z->p->l == z) {
                z->p->l = x;
            } else {
                z->p->r = x;
            }
        }
        if (leftmost == z) {
            if (z->r) {
                leftmost = minimum(x);
            } else {
                leftmost = z->p;
            }
        }
        if (rightmost == z) {
            if (z->l) {
                rightmost = maximum(x);
            } else {
                rightmost = z->p;
            }
        }
    }
    if (y->c != color::red) {
        while ((x != root) && (!x || (x->c == color::black))) {
            if (x == xp->l) {
                base_pointer w = xp->r;
                if (w->c == color::red) {
                    w->c = color::black;
                    xp->c = color::red;
                    rotate_left(xp, root);
                    w = xp->r;
                }
                if ((!w->l || (w->l->c == color::black)) && (!w->r || w->r->c == color::black)) {
                    w->c = color::red;
                    x = xp;
                    xp = xp->p;
                } else {
                    if (!w->r || (w->r->c == color::black)) {
                        w->l->c = color::black;
                        w->c = color::red;
                        rotate_right(w, root);
                        w = xp->r;
                    }
                    w->c = xp->c;
                    xp->c = color::black;
                    if (w->r) {
                        w->r->c = color::black;
                    }
                    rotate_left(xp, root);
                    break;
                }
            } else {
                base_pointer w = xp->l;
                if (w->c == color::red) {
                    w->c = color::black;
                    xp->c = color::red;
                    rotate_right(xp, root);
                    w = xp->l;
                }
                if ((!w->r || (w->r->c == color::black)) && (!w->l || (w->l->c == color::black))) {
                    w->c = color::red;
                    x = xp;
                    xp = xp->p;
                } else {
                    if (!w->l || (w->l->c == color::black)) {
                        w->r->c = color::black;
                        w->c = color::red;
                        rotate_left(w, root);
                        w = xp->l;
                    }
                    w->c = xp->c;
                    xp->c = color::black;
                    if (w->l) {
                        w->l->c = color::black;
                    }
                    rotate_right(xp, root);
                    break;
                }
            }
        }
        if (x) {
            x->c = color::black;
        }
    }
    return y;
}

template< typename value_type, typename compare >
class adapt_compare
{

    static
    const auto & key(const value_type & v)
    {
        return v.k;
    }

    template< typename type >
    static
    const type & key(const type & v)
    {
        return v;
    }

    compare c;

public :

    adapt_compare() = default;

    adapt_compare(const compare & comp)
        : c{comp}
    { ; }

    template< typename L, typename R, typename ...P >
    bool operator () (const L & l, const R & r, P &... p) const
    {
        return c(key(l), key(r), p...);
    }

};

template< typename type >
struct node
        : node_base
{

    union { type value; };

    node() noexcept { ; }

    node(const node &) = delete;
    node(node &&) = delete;
    void operator = (const node &) = delete;
    void operator = (node &&) = delete;

    ~node() { ; }

    type * pointer() noexcept { return &value; }

};

template< typename type >
struct tree_iterator
{

    using value_type = type;
    using reference = type &;
    using pointer = type *;

    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;

    using node_type = node< value_type >;
    using node_pointer = node_type *;

    base_pointer p = nullptr;

    pointer operator -> () const noexcept { return node_pointer(p)->pointer(); }
    reference operator * () const noexcept { return *operator -> (); }

    tree_iterator & operator ++ () noexcept { p = increment(p); return *this; }
    tree_iterator operator ++ (int) noexcept { return {std::exchange(p, increment(p))}; }

    tree_iterator & operator -- () noexcept { p = decrement(p); return *this; }
    tree_iterator operator -- (int) noexcept { return {std::exchange(p, decrement(p))}; }

    bool operator == (const tree_iterator & it) const noexcept { return p == it.p; }
    bool operator != (const tree_iterator & it) const noexcept { return !operator == (it); }

    operator tree_iterator< const type > () const { return {p}; }

};

template< typename type,
          typename compare = std::less< type >,
          typename allocator = std::allocator< type > >
struct tree
{

    using size_type = std::size_t;
    using value_type = type;
    using compare_type = compare;
    using allocator_type = allocator;

private :

    compare_type c;

    using node_type = node< value_type >;
    using node_pointer = node_type *;

    using allocator_traits = typename std::allocator_traits< allocator >::template rebind_traits< node_type >;
    using node_allocator_type = typename allocator_traits::allocator_type;

    node_allocator_type a;

    node_pointer pool = nullptr;

    node_pointer
    get_node()
    {
        if (pool) {
            return std::exchange(pool, node_pointer(pool->r));
        }
        return ::new (allocator_traits::allocate(a, 1)) node_type;
    }

    void put_node(const node_pointer n) noexcept
    {
        n->r = std::exchange(pool, n);
    }

    template< typename ...types >
    node_pointer
    create_node(types &... values)
    {
        const node_pointer n = get_node();
        try {
            allocator_traits::construct(a, n->pointer(), std::forward< types >(values)...);
            return n;
        } catch (...) {
            put_node(n);
            throw;
        }
    }

    void drop_node(const base_pointer n) noexcept
    {
        allocator_traits::destroy(a, node_pointer(n)->pointer());
        put_node(node_pointer(n));
    }

    node_base h{&h};

    size_type s = 0;

    void erase(base_pointer x) noexcept
    {
        while (x) {
            erase(x->r);
            const base_pointer y = x->l;
            drop_node(x);
            x = y;
        }
    }

public :

    tree() = default;

    tree(const tree &) = delete;
    tree(tree &&) = delete;
    void operator = (const tree &) = delete;
    void operator = (tree &&) = delete;

    tree(const compare_type & comp, const allocator_type & alloc)
        : c{comp}
        , a{alloc}
    { ; }

    tree(const compare_type & comp, allocator_type && alloc = allocator_type{})
        : c{comp}
        , a{std::move(alloc)}
    { ; }

    size_type size() const noexcept { return s; }

    bool empty() const noexcept { return (0 == s); }

    void reserve(size_type n)
    {
        node_pointer p = pool;
        while (p && (s < n)) {
            p = p->r;
            --n;
        }
        while (s < n) {
            put_node(get_node());
            --n;
        }
    }

    void shrink_to_fit() noexcept
    {
        while (pool) {
            const base_pointer p = pool->r;
            pool->~node_type();
            allocator_traits::deallocate(a, std::exchange(pool, node_pointer(p)), 1);
        }
    }

    void clear() noexcept
    {
        erase(h.p);
        shrink_to_fit();
        h = {&h};
        s = 0;
    }

    ~tree() noexcept
    {
        clear();
    }

    using iterator = tree_iterator< value_type >;
    using const_iterator = tree_iterator< const value_type >;

    iterator begin() { return {h.l}; }
    iterator end() { return {&h}; }

    const_iterator begin() const { return {h.l}; }
    const_iterator end() const { return {base_pointer(&h)}; }

    const_iterator cbegin() const { return {h.l}; }
    const_iterator cend() const { return {base_pointer(&h)}; }

    iterator
    erase(const const_iterator x) noexcept
    {
        const const_iterator r = std::next(x);
        drop_node(rebalance_for_erase(x.p, h));
        --s;
        return {base_pointer(r.p)};
    }

private :

    static const value_type & value(const base_pointer n) { return *node_pointer(n)->pointer(); }

    static const value_type & value(const iterator x) { return value(x.p); }

    template< typename K, typename ...P >
    pair< base_pointer, base_pointer >
    get_insert_unique_pos(const K & k, P &... p) const
    {
        base_pointer x = h.p;
        auto y = base_pointer(&h);
        bool comp = true;
        while (x) {
            y = x;
            comp = c(k, value(x), p...);
            x = comp ? x->l : x->r;
        }
        base_pointer j = y;
        if (comp) {
            if (j == h.l) {
                return {x, y};
            } else {
                j = decrement(j);
            }
        }
        if (c(value(j), k, p...)) {
            return {x, y};
        }
        return {j, nullptr};
    }

    template< typename K, typename ...P >
    pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(const base_pointer hint, const K & k, P &... p) const
    {
        if (hint == &h) {
            if (!empty() && c(value(h.r), k, p...)) {
                return {nullptr, h.r};
            }
        } else if (c(k, value(hint), p...)) {
            if (hint == h.l) {
                return {hint, hint};
            } else {
                const base_pointer before = decrement(hint);
                if (c(value(before), k, p...)) {
                    if (before->r) {
                        return {hint, hint};
                    } else {
                        return {nullptr, before};
                    }
                }
            }
        } else if (c(value(hint), k, p...)) {
            base_pointer after = hint;
            if (hint == h.r) {
                return {nullptr, hint};
            } else {
                after = increment(after);
                if (c(k, value(after), p...)) {
                    if (hint->r) {
                        return {after, after};
                    } else {
                        return {nullptr, hint};
                    }
                }
            }
        } else {
            return {hint, nullptr};
        }
        return get_insert_unique_pos< K, P... >(k, p...);
    }

    pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(const base_pointer hint) const
    {
        if (hint == &h) {
            if (empty()) {
                return {nullptr, base_pointer(&h)};
            } else {
                return {nullptr, h.r};
            }
        } else {
            if (hint != h.l) {
                const base_pointer before = decrement(hint);
                if (!before->r) {
                    return {nullptr, before};
                }
            }
            return {hint, hint};
        }
    }

    template< typename K, typename ...P >
    base_pointer
    insert_unique(base_pointer l, const base_pointer r, K & k, P &... p)
    {
        if (r) {
            const bool insert_left = (l || (r == &h) || /*c(k, value(r))*/ !c(value(r), k, p...));
            l = create_node< K >(k);
            insert_and_rebalance(insert_left, l, r, h);
            ++s;
            return l;
        }
        return l;
    }

    template< typename K >
    base_pointer
    force_insert_unique(base_pointer l, const base_pointer r, K & k)
    {
        const bool insert_left = (l || (r == &h));
        l = create_node< K >(k);
        insert_and_rebalance(insert_left, l, r, h);
        ++s;
        return l;
    }

public :

    template< typename K = value_type, typename ...P >
    std::pair< iterator, bool >
    insert(K && k, P &... p)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_unique_pos< K, P... >(k, p...);
        return {{insert_unique< K >(lr.k, lr.v, k, p...)}, (lr.v != nullptr)};
    }

    template< typename K = value_type, typename ...P >
    iterator
    insert(const const_iterator hint, K && k, P &... p)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p, k, p...);
        return {insert_unique< K >(lr.k, lr.v, k, p...)};
    }

    template< typename K = value_type >
    iterator
    force_insert(const const_iterator hint, K && k)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p);
        return {force_insert_unique< K >(lr.k, lr.v, k)};
    }

    template< typename K = value_type, typename ...P >
    iterator
    lower_bound(const K & k, P &... p)
    {
        base_pointer l = h.p;
        base_pointer r = &h;
        while (l) {
            if (c(value(l), k, p...)) {
                l = l->r;
            } else {
                r = l;
                l = l->l;
            }
        }
        return {r};
    }

    template< typename K = value_type, typename ...P >
    iterator
    upper_bound(const K & k, P &... p)
    {
        base_pointer l = h.p;
        base_pointer r = &h;
        while (l) {
            if (c(k, value(l), p...)) {
                r = l;
                l = l->l;
            } else {
                l = l->r;
            }
        }
        return {r};
    }

    template< typename K = value_type, typename ...P >
    std::pair< iterator, iterator >
    equal_range(const K & k, P &... p)
    {
        auto l = lower_bound(k, p...);
        auto r = l;
        while ((r != end()) && !c(k, value(r), p...)) {
            ++r;
        }
        return {l, r};
    }

    template< typename K = value_type, typename ...P >
    iterator
    find(const K & k, P &... p)
    {
        const iterator r = lower_bound(k, p...);
        if ((r != end()) && c(k, value(r), p...)) {
            return end();
        }
        return r;
    }

    template< typename C, typename K = value_type, typename ...P >
    iterator
    search(const K & k, C & c, P &... p)
    {
        if (empty()) {
            return end();
        }
        adapt_compare< value_type, C & > a{c};
        base_pointer l = h.p;
        base_pointer r = &h;
        while (l) {
            if (a(k, value(l), p...)) {
                r = l;
                l = l->l;
            } else {
                l = l->r;
            }
        }
        if (r != end()) {
            if (a(k, value(r), p...)) {
                if (r == begin()) {
                    return end();
                }
            } else {
                return r;
            }
        }
        --r;
        if (a(value(r), k, p...)) {
            return end();
        }
        return r;
    }

};

template< typename key_type,
          typename compare = std::less< key_type >,
          typename allocator_type = std::allocator< key_type > >
using set = tree< key_type, compare, allocator_type >;

template< typename key_type,
          typename mapped_type,
          typename compare = std::less< key_type >,
          typename allocator_type = std::allocator< pair< key_type const, mapped_type > > >
using map = tree< typename allocator_type::value_type, adapt_compare< typename allocator_type::value_type, compare >, allocator_type >;

}
