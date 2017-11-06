# sweepline
Fortune's algorithm for Voronoi diagram generating on the plane. Intended for runtime speed and careful handling of corner cases.

How to use:

    using value_type = double;
    struct point 
    { 
        value_type x, y; 
        bool operator < (point const & p) const
        { return std::tie(x, y) < std::tie(p.x, p.y);
    };
    using points = std::vector< point >
    using site = typename points::const_iterator;
    using sweepline_type = sweepline< site, point, value_type >;
    sweepline_type sweepline_{eps};
    points points_;
    // fill points_ with data
    std::sort(std::begin(points_), std::end(points_));
    sweepline_(std::cbegin(points_), std::cend(points_));
    // sweepline_.vertices_ - resulting vertices
    // sweepline_.edges_ - resulting edges
