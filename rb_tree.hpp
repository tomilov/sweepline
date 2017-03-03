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

    node_base(void *)
        : p(nullptr)
        , l(this)
        , r(this)
        , c(color::red)
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
void rotate_left(base_pointer const x, base_pointer & root) noexcept
{
    base_pointer const y = x->r;
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
void rotate_right(base_pointer const x, base_pointer & root) noexcept
{
    base_pointer const y = x->l;
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
void insert_and_rebalance(bool const insert_left,
                          base_pointer x, base_pointer const p,
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
        base_pointer const xpp = x->p->p;
        if (x->p == xpp->l) {
            base_pointer const y = xpp->r;
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
            base_pointer const y = xpp->l;
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
rebalance_for_erase(base_pointer const z, node_base & h) noexcept
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

template< typename type >
struct node;

template< typename type >
struct node< type const >
        : node_base
{

    type const * pointer() const noexcept { return &storage.value; }

protected :

    union storage_type
    {

        storage_type() noexcept { ; }
        ~storage_type() noexcept { ; }

        type value;

    } storage;

};

template< typename type >
struct node
        : node< type const >
{

    type * pointer() noexcept { return &node< type const >::storage.value; }

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

    bool operator == (tree_iterator const & it) const noexcept { return p == it.p; }
    bool operator != (tree_iterator const & it) const noexcept { return !operator == (it); }

    operator tree_iterator< type const > () const
    {
        return {p};
    }

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

    node_pointer get_node()
    {
        if (pool) {
            return std::exchange(pool, node_pointer(pool->r));
        }
        return ::new (allocator_traits::allocate(a, 1)) node_type;
    }

    void put_node(node_pointer const n) noexcept
    {
        n->r = std::exchange(pool, n);
    }

    template< typename... types >
    node_pointer
    create_node(types &&... values)
    {
        node_pointer const n = get_node();
        try {
            allocator_traits::construct(a, n->pointer(), std::forward< types >(values)...);
            return n;
        } catch (...) {
            put_node(n);
            throw;
        }
    }

    void drop_node(base_pointer const n) noexcept
    {
        allocator_traits::destroy(a, node_pointer(n)->pointer());
        put_node(node_pointer(n));
    }

    node_base h{{}};

    size_type s = 0;

    void erase(base_pointer x) noexcept
    {
        while (x) {
            erase(x->r);
            base_pointer const y = x->l;
            drop_node(x);
            x = y;
        }
    }

public :

    tree() = default;

    tree(compare_type const & comp, allocator_type const & alloc)
        : c(comp)
        , a(alloc)
    { ; }

    explicit
    tree(compare_type const & comp, allocator_type && alloc = allocator_type{})
        : c(comp)
        , a(std::move(alloc))
    { ; }

    size_type size() const noexcept { return s; }

    bool empty() const noexcept { return (0 == s); }

    void reserve(size_type n)
    {
        node_pointer p = pool;
        while (s < n) {
            if (p) {
                p = p->r;
            } else {
                while (s < n) {
                    put_node(get_node());
                    --n;
                }
                return;
            }
            --n;
        }
    }

    void shrink_to_fit() noexcept
    {
        while (pool) {
            base_pointer const p = pool->r;
            pool->~node_type();
            allocator_traits::deallocate(a, std::exchange(pool, node_pointer(p)), 1);
        }
    }

    void clear() noexcept
    {
        erase(h.p);
        h.c = color::red;
        h.p = nullptr;
        h.l = &h;
        h.r = &h;
        s = 0;
        shrink_to_fit();
    }

    ~tree() noexcept
    {
        clear();
    }

    tree(tree const &) = delete;
    tree(tree &&) = delete;
    void operator = (tree const &) = delete;
    void operator = (tree &&) = delete;

    using iterator = tree_iterator< value_type >;
    using const_iterator = tree_iterator< value_type const >;

    iterator begin() { return {h.l}; }
    iterator end() { return {base_pointer(&h)}; }

    const_iterator begin() const { return {h.l}; }
    const_iterator end() const { return {base_pointer(&h)}; }

    const_iterator cbegin() const { return {h.l}; }
    const_iterator cend() const { return {base_pointer(&h)}; }

    iterator
    erase(iterator const x) noexcept
    {
        iterator const r = std::next(x);
        drop_node(rebalance_for_erase(x.p, h));
        --s;
        return r;
    }

private :

    static value_type const & value(base_pointer const n) { return *node_pointer(n)->pointer(); }

    pair< base_pointer, base_pointer >
    get_insert_unique_pos(value_type const & v) const
    {
        base_pointer x = h.p;
        auto y = base_pointer(&h);
        bool comp = true;
        while (x) {
            y = x;
            comp = c(v, value(x));
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
        if (c(value(j), v)) {
            return {x, y};
        }
        return {j, nullptr};
    }

    template< typename K >
    pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(base_pointer const hint, K const & k) const
    {
        if (hint == &h) {
            if (!empty() && c(value(h.r), k)) {
                return {nullptr, h.r};
            } else {
                return get_insert_unique_pos(k);
            }
        } else if (c(k, value(hint))) {
            if (hint == h.l) {
                return {hint, hint};
            } else {
                base_pointer const before = decrement(hint);
                if (c(value(before), k)) {
                    if (before->r) {
                        return {hint, hint};
                    } else {
                        return {nullptr, before};
                    }
                } else {
                    return get_insert_unique_pos(k);
                }
            }
        } else if (c(value(hint), k)) {
            base_pointer after = hint;
            if (hint == h.r) {
                return {nullptr, hint};
            } else {
                after = increment(after);
                if (c(k, value(after))) {
                    if (hint->r) {
                        return {after, after};
                    } else {
                        return {nullptr, hint};
                    }
                } else {
                    return get_insert_unique_pos(k);
                }
            }
        } else {
            return {hint, nullptr};
        }
    }

    pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(base_pointer const hint) const
    {
        if (hint == &h) {
            if (empty()) {
                return {nullptr, base_pointer(&h)};
            } else {
                return {nullptr, h.r};
            }
        } else {
            if (hint == h.l) {
                return {hint, hint};
            } else {
                base_pointer const before = decrement(hint);
                if (before->r) {
                    return {hint, hint};
                } else {
                    return {nullptr, before};
                }
            }
        }
    }

    template< typename K >
    base_pointer
    insert_unique(base_pointer const l, base_pointer const r, K && k)
    {
        if (r) {
            bool const insert_left = (l || (r == &h) || /*c(k, value(r))*/ !c(value(r), k));
            node_pointer const z = create_node(std::forward< K >(k));
            insert_and_rebalance(insert_left, z, r, h);
            ++s;
            return z;
        }
        return l;
    }

    template< typename K >
    base_pointer
    force_insert_unique(base_pointer const l, base_pointer const r, K && k)
    {
        bool const insert_left = (l || (r == &h));
        node_pointer const n = create_node(std::forward< K >(k));
        insert_and_rebalance(insert_left, n, r, h);
        ++s;
        return n;
    }

public :

    template< typename K = value_type >
    std::pair< iterator, bool >
    insert(K && k)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_unique_pos(k);
        return {{insert_unique(lr.k, lr.v, std::forward< K >(k))}, (lr.v != nullptr)};
    }

    template< typename K = value_type >
    iterator
    insert(iterator const hint, K && k)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p, k);
        return {insert_unique(lr.k, lr.v, std::forward< K >(k))};
    }

    template< typename K = value_type >
    iterator
    force_insert(iterator const hint, K && k)
    {
        pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p);
        return {force_insert_unique(lr.k, lr.v, std::forward< K >(k))};
    }

    template< typename K = value_type >
    iterator
    lower_bound(K const & k)
    {
        base_pointer l = h.p;
        base_pointer r = &h;
        while (l) {
            if (c(value(l), k)) {
                l = l->r;
            } else {
                r = l;
                l = l->l;
            }
        }
        return {r};
    }

    template< typename K = value_type >
    iterator
    find(K const & k)
    {
        iterator const r = lower_bound(k);
        if ((r != end()) && c(k, *r)) {
            return end();
        }
        return r;
    }

};

template< typename key_type,
          typename compare = std::less< key_type >,
          typename allocator_type = std::allocator< key_type > >
using set = tree< key_type, compare, allocator_type >;

template< typename value_type, typename compare >
class adapt_compare
{

    static
    auto const & key(value_type const & v)
    {
        return v.k;
    }

    template< typename type >
    static
    type const & key(type const & v)
    {
        return v;
    }

    compare c;

public :

    adapt_compare(compare const & comp)
        : c(comp)
    { ; }

    template< typename L, typename R >
    bool operator () (L const & l, R const & r) const
    {
        return c(key(l), key(r));
    }

};

template< typename key_type,
          typename mapped_type,
          typename compare = std::less< key_type >,
          typename allocator_type = std::allocator< pair< key_type const, mapped_type > > >
using map = tree< typename allocator_type::value_type, adapt_compare< typename allocator_type::value_type, compare >, allocator_type >;

}
