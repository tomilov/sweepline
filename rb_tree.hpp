#pragma once

#include <type_traits>
#include <utility>
#include <iterator>
#include <memory>

#include <cassert>

namespace rb_tree
{

enum class color : bool { red = false, black = true };

struct node_base;

using base_pointer = node_base *;

struct node_base
{

    base_pointer p = nullptr;
    base_pointer l = this;
    base_pointer r = this;
    color c = color::red;

#ifdef DEBUG
    ~node_base()
    {
        p = nullptr;
        l = nullptr;
        r = nullptr;
    }
#endif

};

inline base_pointer minimum(base_pointer x) noexcept { while (x->l) x = x->l; return x; }
inline base_pointer maximum(base_pointer x) noexcept { while (x->r) x = x->r; return x; }

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
        base_pointer y = x->l;
        while (y->r) {
            y = y->r;
        }
        x = y;
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

template< typename type >
struct node : node_base
{

    using value_type = type;

    std::aligned_storage_t< sizeof(type), alignof(type) > storage;

    type * pointer() noexcept { return reinterpret_cast< type * >(static_cast< void * >(&storage)); }
    type const * pointer() const noexcept { return reinterpret_cast< type const * >(static_cast< void const * >(&storage)); }

};

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
                          base_pointer x,
                          base_pointer const p,
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
rebalance_for_erase(base_pointer const z,
                    node_base & h) noexcept
{
    base_pointer & root = h.p;
    base_pointer & leftmost = h.l;
    base_pointer & rightmost = h.r;
    base_pointer y = z;
    base_pointer x = nullptr;
    base_pointer xp = nullptr;
    if (y->l) {
        if (y->r) {
            y = y->r;
            while (y->l) {
                y = y->l;
            }
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

    void put_node(node_pointer const p) noexcept
    {
        // p->~node_type(); // on pool_free()
        p->r = std::exchange(pool, p);
    }

    void pool_free() noexcept
    {
        while (pool) {
            base_pointer const p = pool->r;
            pool->~node_type();
            allocator_traits::deallocate(a, std::exchange(pool, node_pointer(p)), 1);
        }
    }

    template< typename ...Args >
    void construct_node(node_pointer const n, Args && ...args)
    {
        try {
            allocator_traits::construct(a, n->pointer(), std::forward< Args >(args)...);
        } catch (...) {
            put_node(n);
            throw;
        }
    }

    template< typename... Args >
    node_pointer
    create_node(Args &&... args)
    {
        node_pointer const t = get_node();
        construct_node(t, std::forward< Args >(args)...);
        return t;
    }

    void destroy_node(node_pointer const p) noexcept
    {
        allocator_traits::destroy(a, p->pointer());
    }

    void drop_node(node_pointer const p) noexcept
    {
        destroy_node(p);
        put_node(p);
    }

    node_base h;

    size_type s = 0;

    void erase(base_pointer x)
    {
        while (x) {
            erase(x->r);
            base_pointer const y = x->l;
            drop_node(node_pointer(x));
            x = y;
        }
    }

public :

    tree() = default;

    tree(compare_type const & _c,
         allocator_type const & _a)
        : c(_c)
        , a(node_allocator_type(_a))
    { ; }

    tree(compare_type const & _c,
         allocator_type && _a = allocator_type())
        : c(_c)
        , a(node_allocator_type(std::move(_a)))
    { ; }

    size_type size() const { return s; }

    bool empty() const { return (0 == size()); }

    void reserve(size_type n)
    {
        while (size() < n) {
            put_node(get_node());
            --n;
        }
    }

    void clear()
    {
        erase(h.p);
        h.c = color::red;
        h.p = nullptr;
        h.l = &h;
        h.r = &h;
        s = 0;
        pool_free();
    }

    ~tree() noexcept
    {
        clear();
    }

    using iterator = tree_iterator< value_type >;

    iterator begin() const { return {h.l}; }
    iterator end() const { return {base_pointer(&h)}; }

    iterator
    erase(iterator const p)
    {
        iterator const r = std::next(p);
        drop_node(node_pointer(rebalance_for_erase(p.p, h)));
        --s;
        return r;
    }

private :

    static value_type const & value(base_pointer const p) { return *node_pointer(p)->pointer(); }

    std::pair< base_pointer, base_pointer >
    get_insert_unique_pos(value_type const & v)
    {
        base_pointer x = h.p;
        base_pointer y = &h;
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
    std::pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(base_pointer const hint, K const & v)
    {
        if (hint == &h) {
            if (!empty() && c(value(h.r), v)) {
                return {nullptr, h.r};
            } else {
                return get_insert_unique_pos(v);
            }
        } else if (c(v, value(hint))) {
            if (hint == h.l) {
                return {hint, hint};
            } else {
                base_pointer const before = decrement(hint);
                if (c(value(before), v)) {
                    if (before->r) {
                        return {hint, hint};
                    } else {
                        return {nullptr, before};
                    }
                } else {
                    return get_insert_unique_pos(v);
                }
            }
        } else if (c(value(hint), v)) {
            base_pointer after = hint;
            if (hint == h.r) {
                return {nullptr, hint};
            } else {
                after = increment(after);
                if (c(v, value(after))) {
                    if (hint->r) {
                        return {after, after};
                    } else {
                        return {nullptr, hint};
                    }
                } else {
                    return get_insert_unique_pos(v);
                }
            }
        } else {
            return {hint, nullptr};
        }
    }

    template< typename K >
    base_pointer
    insert_unique(base_pointer const l, base_pointer const r, K && v)
    {
        if (r) {
            bool const insert_left = (l || (r == &h) || /*c(v, value(r))*/ !c(value(r), v));
            node_pointer const z = create_node(std::forward< K >(v));
            insert_and_rebalance(insert_left, z, r, h);
            ++s;
            return z;
        }
        return l;
    }

    std::pair< base_pointer, base_pointer >
    get_insert_hint_unique_pos(base_pointer const hint)
    {
        if (hint == &h) {
            if (empty()) {
                return {nullptr, &h};
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

    base_pointer
    insert_unique(base_pointer const l, base_pointer const r, node_pointer const z)
    {
        bool const insert_left = (l || (r == &h));
        insert_and_rebalance(insert_left, z, r, h);
        ++s;
        return z;
    }

public :

    template< typename K >
    iterator
    insert(K && v)
    {
        std::pair< base_pointer, base_pointer > const lr = get_insert_unique_pos(v);
        return {insert_unique(lr.first, lr.second, std::forward< K >(v))};
    }

    template< typename K = value_type >
    iterator
    insert(iterator const hint, K && v)
    {
        std::pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p, v);
        return {insert_unique(lr.first, lr.second, std::forward< K >(v))};
    }

    template< typename K = value_type >
    iterator
    force_insert(iterator const hint, K && v)
    {
        std::pair< base_pointer, base_pointer > const lr = get_insert_hint_unique_pos(hint.p);
        return {insert_unique(lr.first, lr.second, create_node(std::forward< K >(v)))};
    }

    template< typename K >
    iterator
    lower_bound(K && k)
    {
        base_pointer x = h.p;
        base_pointer y = &h;
        while (x) {
            if (c(value(x), k)) {
                x = x->r;
            } else {
                y = x;
                x = x->l;
            }
        }
        return {y};
    }

};

template< typename first_type, typename second_type >
struct pair
{

    first_type first;
    second_type second;

    operator first_type const & () const
    {
        return first;
    }

};

template< typename key_type, typename value_type, typename key_compare_type = std::less< key_type > >
using map = tree< pair< key_type, value_type >, key_compare_type >;

}
