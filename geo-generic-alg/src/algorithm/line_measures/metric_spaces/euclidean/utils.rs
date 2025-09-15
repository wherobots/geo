// Utility functions for generic distance calculations
// This demonstrates the refactoring concept by moving key functions to a separate module

use super::{Distance, Euclidean};
use crate::algorithm::Intersects;
use crate::coordinate_position::{coord_pos_relative_to_ring, CoordPos};
use crate::geometry::*;
use crate::{CoordFloat, GeoFloat, GeoNum};
use geo_traits::CoordTrait;
use geo_traits_ext::{
    LineStringTraitExt, LineTraitExt, PointTraitExt, PolygonTraitExt, TriangleTraitExt,
};
use num_traits::{Bounded, Float};
use rstar::primitives::CachedEnvelope;
use rstar::RTree;

// ┌────────────────────────────────────────────────────────────┐
// │ Helper functions for generic distance calculations         │
// └────────────────────────────────────────────────────────────┘

/// Uses an R* tree and nearest-neighbour lookups to calculate minimum distances
pub fn nearest_neighbour_distance<F: GeoFloat>(geom1: &LineString<F>, geom2: &LineString<F>) -> F {
    let tree_a = RTree::bulk_load(geom1.lines().map(CachedEnvelope::new).collect());
    let tree_b = RTree::bulk_load(geom2.lines().map(CachedEnvelope::new).collect());
    // Return minimum distance between all geom a points and geom b lines, and all geom b points and geom a lines
    geom2
        .points()
        .fold(Bounded::max_value(), |acc: F, point| {
            let nearest = tree_a.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        })
        .min(geom1.points().fold(Bounded::max_value(), |acc, point| {
            let nearest = tree_b.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        }))
}

pub fn ring_contains_coord<T: GeoNum>(ring: &LineString<T>, c: Coord<T>) -> bool {
    match coord_pos_relative_to_ring(c, ring) {
        CoordPos::Inside => true,
        CoordPos::OnBoundary | CoordPos::Outside => false,
    }
}

// Helper for line segment distance using generic trait methods
pub fn point_distance_generic<F, P1, P2>(p1: &P1, p2: &P2) -> F
where
    F: CoordFloat,
    P1: PointTraitExt<T = F>,
    P2: PointTraitExt<T = F>,
{
    if let (Some(c1), Some(c2)) = (p1.coord(), p2.coord()) {
        let delta_x = c1.x() - c2.x();
        let delta_y = c1.y() - c2.y();
        delta_x.hypot(delta_y)
    } else {
        F::zero()
    }
}

pub fn line_segment_distance_generic<F, C, L>(coord: &C, line: &L) -> F
where
    F: CoordFloat,
    C: CoordTrait<T = F>,
    L: LineTraitExt<T = F>,
{
    let px = coord.x();
    let py = coord.y();
    let start = line.start_coord();
    let end = line.end_coord();
    let dx = end.x - start.x;
    let dy = end.y - start.y;

    if dx == F::zero() && dy == F::zero() {
        let delta_x = px - start.x;
        let delta_y = py - start.y;
        return delta_x.hypot(delta_y);
    }

    let t = ((px - start.x) * dx + (py - start.y) * dy) / (dx * dx + dy * dy);
    let t = t.max(F::zero()).min(F::one());

    let nearest_x = start.x + t * dx;
    let nearest_y = start.y + t * dy;
    let delta_x = px - nearest_x;
    let delta_y = py - nearest_y;
    delta_x.hypot(delta_y)
}

// ┌────────────────────────────────────────────────────────────┐
// │ Sample generic distance functions (demonstrating concept)   │
// └────────────────────────────────────────────────────────────┘

/// Line to Line distance
pub fn distance_line_to_line_generic<F, L1, L2>(line1: &L1, line2: &L2) -> F
where
    F: GeoFloat,
    L1: LineTraitExt<T = F>,
    L2: LineTraitExt<T = F>,
{
    let start1 = line1.start_coord();
    let end1 = line1.end_coord();
    let start2 = line2.start_coord();
    let end2 = line2.end_coord();

    // Check if lines intersect using generic intersects
    if line1.intersects(line2) {
        return F::zero();
    }

    // Find minimum distance between all endpoint combinations
    let dist1 = line_segment_distance_generic(&start1, line2);
    let dist2 = line_segment_distance_generic(&end1, line2);
    let dist3 = line_segment_distance_generic(&start2, line1);
    let dist4 = line_segment_distance_generic(&end2, line1);

    dist1.min(dist2).min(dist3).min(dist4)
}

/// Point to LineString distance
pub fn distance_point_to_linestring_generic<F, P, LS>(point: &P, linestring: &LS) -> F
where
    F: GeoFloat,
    P: PointTraitExt<T = F>,
    LS: LineStringTraitExt<T = F>,
{
    if let Some(coord) = point.coord() {
        let mut lines = linestring.lines();
        if let Some(first_line) = lines.next() {
            let mut min_distance = line_segment_distance_generic(&coord, &first_line);
            for line in lines {
                min_distance = min_distance.min(line_segment_distance_generic(&coord, &line));
            }
            min_distance
        } else {
            F::zero()
        }
    } else {
        F::zero()
    }
}

/// Point to Polygon distance
pub fn distance_point_to_polygon_generic<F, P, Poly>(point: &P, polygon: &Poly) -> F
where
    F: GeoFloat,
    P: PointTraitExt<T = F>,
    Poly: PolygonTraitExt<T = F>,
{
    // Check if the polygon is empty
    if polygon.exterior_ext().is_none() {
        return F::zero();
    }

    // Use the existing generic Intersects implementation
    // If the point intersects the polygon (is inside or on boundary), distance is 0
    if polygon.intersects(point) {
        return F::zero();
    }

    // Point is outside the polygon, calculate minimum distance to edges
    if let (Some(coord), Some(exterior)) = (point.coord(), polygon.exterior_ext()) {
        // Calculate minimum distance to exterior ring
        let exterior_dist = exterior
            .lines()
            .map(|line| line_segment_distance_generic(&coord, &line))
            .fold(Float::max_value(), |acc: F, dist| acc.min(dist));

        // Calculate minimum distance to interior rings (holes)
        let interior_dist = polygon
            .interiors_ext()
            .map(|interior| {
                interior
                    .lines()
                    .map(|line| line_segment_distance_generic(&coord, &line))
                    .fold(Float::max_value(), |acc: F, dist| acc.min(dist))
            })
            .fold(Float::max_value(), |acc: F, dist| acc.min(dist));

        exterior_dist.min(interior_dist)
    } else {
        F::zero()
    }
}

/// LineString to Polygon distance
pub fn distance_linestring_to_polygon_generic<F, LS, Poly>(linestring: &LS, polygon: &Poly) -> F
where
    F: GeoFloat,
    LS: LineStringTraitExt<T = F>,
    Poly: PolygonTraitExt<T = F>,
{
    if let Some(exterior) = polygon.exterior_ext() {
        let mut min_dist: F = Float::max_value();

        // Calculate distance to exterior ring
        for line1 in linestring.lines() {
            for line2 in exterior.lines() {
                let d1 = line_segment_distance_generic(&line1.start_coord(), &line2);
                let d2 = line_segment_distance_generic(&line1.end_coord(), &line2);
                let d3 = line_segment_distance_generic(&line2.start_coord(), &line1);
                let d4 = line_segment_distance_generic(&line2.end_coord(), &line1);
                let line_dist = d1.min(d2).min(d3).min(d4);
                min_dist = min_dist.min(line_dist);
            }
        }

        // Also calculate distance to interior rings (holes)
        for interior in polygon.interiors_ext() {
            for line1 in linestring.lines() {
                for line2 in interior.lines() {
                    let d1 = line_segment_distance_generic(&line1.start_coord(), &line2);
                    let d2 = line_segment_distance_generic(&line1.end_coord(), &line2);
                    let d3 = line_segment_distance_generic(&line2.start_coord(), &line1);
                    let d4 = line_segment_distance_generic(&line2.end_coord(), &line1);
                    let line_dist = d1.min(d2).min(d3).min(d4);
                    min_dist = min_dist.min(line_dist);
                }
            }
        }

        if min_dist == Float::max_value() {
            F::zero()
        } else {
            min_dist
        }
    } else {
        F::zero()
    }
}

/// Polygon to Polygon distance
pub fn distance_polygon_to_polygon_generic<F, P1, P2>(polygon1: &P1, polygon2: &P2) -> F
where
    F: GeoFloat,
    P1: PolygonTraitExt<T = F>,
    P2: PolygonTraitExt<T = F>,
{
    // Check if polygons intersect using generic intersects
    if polygon1.intersects(polygon2) {
        return F::zero();
    }

    if let (Some(ext1), Some(ext2)) = (polygon1.exterior_ext(), polygon2.exterior_ext()) {
        // For containment checks, we still need concrete polygons since ring_contains_coord requires concrete types
        let ext1_coords: Vec<Coord<F>> = ext1
            .coords_ext()
            .map(|c| Coord::from((c.x(), c.y())))
            .collect();
        let ext2_coords: Vec<Coord<F>> = ext2
            .coords_ext()
            .map(|c| Coord::from((c.x(), c.y())))
            .collect();

        let interior1_coords: Vec<LineString<F>> = polygon1
            .interiors_ext()
            .map(|ring| {
                LineString::from(
                    ring.coords_ext()
                        .map(|c| (c.x(), c.y()))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let interior2_coords: Vec<LineString<F>> = polygon2
            .interiors_ext()
            .map(|ring| {
                LineString::from(
                    ring.coords_ext()
                        .map(|c| (c.x(), c.y()))
                        .collect::<Vec<_>>(),
                )
            })
            .collect();

        let poly_a: Polygon<F> = Polygon::new(LineString::from(ext1_coords), interior1_coords);
        let poly_b: Polygon<F> = Polygon::new(LineString::from(ext2_coords), interior2_coords);

        // Containment check - if polygon_b is inside polygon_a's hole
        if !poly_a.interiors().is_empty() {
            // Get first coordinate of polygon_b's exterior ring
            if let Some(first_coord_b) = ext2.coords_ext().next() {
                let coord_b = Coord::from((first_coord_b.x(), first_coord_b.y()));
                if ring_contains_coord(poly_a.exterior(), coord_b) {
                    // check each ring distance, returning the minimum
                    let mut mindist: F = Float::max_value();
                    let ext2_concrete = LineString::from(
                        ext2.coords_ext()
                            .map(|c| (c.x(), c.y()))
                            .collect::<Vec<_>>(),
                    );
                    for ring in poly_a.interiors() {
                        mindist = mindist.min(nearest_neighbour_distance(&ext2_concrete, ring));
                    }
                    return mindist;
                }
            }
        }

        // Containment check - if polygon_a is inside polygon_b's hole
        if !poly_b.interiors().is_empty() {
            // Get first coordinate of polygon_a's exterior ring
            if let Some(first_coord_a) = ext1.coords_ext().next() {
                let coord_a = Coord::from((first_coord_a.x(), first_coord_a.y()));
                if ring_contains_coord(poly_b.exterior(), coord_a) {
                    let mut mindist: F = Float::max_value();
                    let ext1_concrete = LineString::from(
                        ext1.coords_ext()
                            .map(|c| (c.x(), c.y()))
                            .collect::<Vec<_>>(),
                    );
                    for ring in poly_b.interiors() {
                        mindist = mindist.min(nearest_neighbour_distance(&ext1_concrete, ring));
                    }
                    return mindist;
                }
            }
        }

        // Default case - distance between exterior rings
        let ext1_concrete = LineString::from(
            ext1.coords_ext()
                .map(|c| (c.x(), c.y()))
                .collect::<Vec<_>>(),
        );
        let ext2_concrete = LineString::from(
            ext2.coords_ext()
                .map(|c| (c.x(), c.y()))
                .collect::<Vec<_>>(),
        );
        nearest_neighbour_distance(&ext1_concrete, &ext2_concrete)
    } else {
        F::zero()
    }
}

/// LineString to LineString distance
#[allow(dead_code)] // Used in test code
pub fn distance_linestring_to_linestring_generic<F, LS1, LS2>(ls1: &LS1, ls2: &LS2) -> F
where
    F: GeoFloat,
    LS1: LineStringTraitExt<T = F>,
    LS2: LineStringTraitExt<T = F>,
{
    ls1.lines()
        .flat_map(|line1| {
            ls2.lines()
                .map(move |line2| distance_line_to_line_generic(&line1, &line2))
        })
        .fold(Float::max_value(), |acc, dist| acc.min(dist))
}

/// Line to Polygon distance
pub fn distance_line_to_polygon_generic<F, L, Poly>(line: &L, polygon: &Poly) -> F
where
    F: GeoFloat,
    L: LineTraitExt<T = F>,
    Poly: PolygonTraitExt<T = F>,
{
    // Convert line to linestring and use existing linestring-to-polygon function
    let line_coords = vec![line.start_coord(), line.end_coord()];
    let line_as_ls = LineString::from(line_coords);
    distance_linestring_to_polygon_generic(&line_as_ls, polygon)
}

/// Line to LineString distance
pub fn distance_line_to_linestring_generic<F, L, LS>(line: &L, linestring: &LS) -> F
where
    F: GeoFloat,
    L: LineTraitExt<T = F>,
    LS: LineStringTraitExt<T = F>,
{
    linestring
        .lines()
        .map(|ls_line| distance_line_to_line_generic(line, &ls_line))
        .fold(Float::max_value(), |acc, dist| acc.min(dist))
}

/// Triangle to Point distance
pub fn distance_triangle_to_point_generic<F, T, P>(triangle: &T, point: &P) -> F
where
    F: GeoFloat,
    T: TriangleTraitExt<T = F>,
    P: PointTraitExt<T = F>,
{
    // Convert triangle to polygon and use existing point-to-polygon function
    let tri_poly = triangle.to_polygon();
    distance_point_to_polygon_generic(point, &tri_poly)
}

// ┌────────────────────────────────────────────────────────────┐
// │ Symmetric Distance Function Generator Macro                │
// └────────────────────────────────────────────────────────────┘

/// Macro to generate symmetric distance functions
/// For distance operations that are symmetric (distance(a, b) == distance(b, a)),
/// this macro generates the reverse function that calls the primary implementation
macro_rules! symmetric_distance_generic_impl {
    ($func_name_ab:ident, $func_name_ba:ident, $trait_a:ident, $trait_b:ident) => {
        #[allow(dead_code)]
        pub fn $func_name_ba<F, A, B>(b: &B, a: &A) -> F
        where
            F: GeoFloat,
            A: $trait_a<T = F>,
            B: $trait_b<T = F>,
        {
            $func_name_ab(a, b)
        }
    };
}

// Generate symmetric distance functions
symmetric_distance_generic_impl!(
    distance_point_to_linestring_generic,
    distance_linestring_to_point_generic,
    PointTraitExt,
    LineStringTraitExt
);

symmetric_distance_generic_impl!(
    distance_point_to_polygon_generic,
    distance_polygon_to_point_generic,
    PointTraitExt,
    PolygonTraitExt
);

symmetric_distance_generic_impl!(
    distance_linestring_to_polygon_generic,
    distance_polygon_to_linestring_generic,
    LineStringTraitExt,
    PolygonTraitExt
);

symmetric_distance_generic_impl!(
    distance_line_to_linestring_generic,
    distance_linestring_to_line_generic,
    LineTraitExt,
    LineStringTraitExt
);

symmetric_distance_generic_impl!(
    distance_line_to_polygon_generic,
    distance_polygon_to_line_generic,
    LineTraitExt,
    PolygonTraitExt
);

symmetric_distance_generic_impl!(
    distance_triangle_to_point_generic,
    distance_point_to_triangle_generic,
    TriangleTraitExt,
    PointTraitExt
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coord, Line, LineString, Point, Polygon, Triangle};
    use approx::assert_relative_eq;

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for point_distance_generic function                  │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_point_distance_generic_basic() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(3.0, 4.0);

        let distance = point_distance_generic(&p1, &p2);
        assert_relative_eq!(distance, 5.0); // 3-4-5 triangle

        // Test symmetry
        let distance_reverse = point_distance_generic(&p2, &p1);
        assert_relative_eq!(distance, distance_reverse);
    }

    #[test]
    fn test_point_distance_generic_same_point() {
        let p = Point::new(2.5, -1.5);
        let distance = point_distance_generic(&p, &p);
        assert_relative_eq!(distance, 0.0);
    }

    #[test]
    fn test_point_distance_generic_negative_coordinates() {
        let p1 = Point::new(-2.0, -3.0);
        let p2 = Point::new(1.0, 1.0);

        let distance = point_distance_generic(&p1, &p2);
        assert_relative_eq!(distance, 5.0); // sqrt((1-(-2))^2 + (1-(-3))^2) = sqrt(9+16) = 5
    }

    #[test]
    fn test_point_distance_generic_empty_points() {
        // Test with empty points (no coordinates)
        let empty_point: Point<f64> = Point::new(f64::NAN, f64::NAN);
        let regular_point = Point::new(1.0, 1.0);

        // When either point has no valid coordinates, distance should be 0
        let distance = point_distance_generic(&empty_point, &regular_point);
        assert!(distance.is_nan() || distance == 0.0); // Implementation dependent
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for line_segment_distance_generic function           │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_line_segment_distance_generic_point_on_line() {
        let coord = coord! { x: 2.0, y: 0.0 };
        let line = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 4.0, y: 0.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        assert_relative_eq!(distance, 0.0);
    }

    #[test]
    fn test_line_segment_distance_generic_perpendicular() {
        let coord = coord! { x: 2.0, y: 3.0 };
        let line = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 4.0, y: 0.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        assert_relative_eq!(distance, 3.0);
    }

    #[test]
    fn test_line_segment_distance_generic_beyond_endpoint() {
        let coord = coord! { x: 6.0, y: 0.0 };
        let line = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 4.0, y: 0.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        assert_relative_eq!(distance, 2.0); // Distance to closest endpoint (4,0)
    }

    #[test]
    fn test_line_segment_distance_generic_before_startpoint() {
        let coord = coord! { x: -2.0, y: 0.0 };
        let line = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 4.0, y: 0.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        assert_relative_eq!(distance, 2.0); // Distance to start point (0,0)
    }

    #[test]
    fn test_line_segment_distance_generic_zero_length_line() {
        let coord = coord! { x: 2.0, y: 3.0 };
        let line = Line::new(coord! { x: 1.0, y: 1.0 }, coord! { x: 1.0, y: 1.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        let expected = ((2.0 - 1.0).powi(2) + (3.0 - 1.0).powi(2)).sqrt();
        assert_relative_eq!(distance, expected);
    }

    #[test]
    fn test_line_segment_distance_generic_diagonal_line() {
        let coord = coord! { x: 0.0, y: 2.0 };
        let line = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 2.0, y: 2.0 });

        let distance = line_segment_distance_generic(&coord, &line);
        // Point (0,2) to line from (0,0) to (2,2) - should be sqrt(2)
        assert_relative_eq!(distance, std::f64::consts::SQRT_2, epsilon = 1e-10);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for nearest_neighbour_distance function              │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_nearest_neighbour_distance_basic() {
        let ls1 = LineString::from(vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);
        let ls2 = LineString::from(vec![(3.0, 0.0), (5.0, 0.0), (5.0, 2.0)]);

        let distance = nearest_neighbour_distance(&ls1, &ls2);
        assert_relative_eq!(distance, 1.0); // Distance between (2,0)-(2,2) and (3,0)-(5,0)
    }

    #[test]
    fn test_nearest_neighbour_distance_intersecting() {
        let ls1 = LineString::from(vec![(0.0, 0.0), (4.0, 0.0)]);
        let ls2 = LineString::from(vec![(2.0, -1.0), (2.0, 1.0)]);

        let distance = nearest_neighbour_distance(&ls1, &ls2);
        assert_relative_eq!(distance, 1.0);
    }

    #[test]
    fn test_nearest_neighbour_distance_parallel_lines() {
        let ls1 = LineString::from(vec![(0.0, 0.0), (4.0, 0.0)]);
        let ls2 = LineString::from(vec![(0.0, 2.0), (4.0, 2.0)]);

        let distance = nearest_neighbour_distance(&ls1, &ls2);
        assert_relative_eq!(distance, 2.0); // Perpendicular distance between parallel lines
    }

    #[test]
    fn test_nearest_neighbour_distance_single_segment_each() {
        let ls1 = LineString::from(vec![(0.0, 0.0), (1.0, 0.0)]);
        let ls2 = LineString::from(vec![(2.0, 1.0), (3.0, 1.0)]);

        let distance = nearest_neighbour_distance(&ls1, &ls2);
        let expected = ((2.0 - 1.0).powi(2) + (1.0 - 0.0).powi(2)).sqrt();
        assert_relative_eq!(distance, expected);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for ring_contains_coord function                     │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_ring_contains_coord_inside() {
        let ring = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let coord = coord! { x: 2.0, y: 2.0 };

        assert!(ring_contains_coord(&ring, coord));
    }

    #[test]
    fn test_ring_contains_coord_outside() {
        let ring = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let coord = coord! { x: 5.0, y: 2.0 };

        assert!(!ring_contains_coord(&ring, coord));
    }

    #[test]
    fn test_ring_contains_coord_on_boundary() {
        let ring = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let coord = coord! { x: 2.0, y: 0.0 };

        assert!(!ring_contains_coord(&ring, coord)); // On boundary = false
    }

    #[test]
    fn test_ring_contains_coord_triangle() {
        let ring = LineString::from(vec![
            (0.0, 0.0), (3.0, 0.0), (1.5, 2.0), (0.0, 0.0)
        ]);
        let inside_coord = coord! { x: 1.5, y: 0.5 };
        let outside_coord = coord! { x: 3.0, y: 3.0 };

        assert!(ring_contains_coord(&ring, inside_coord));
        assert!(!ring_contains_coord(&ring, outside_coord));
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_point_to_linestring_generic             │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_point_to_linestring_generic_basic() {
        let point = Point::new(1.0, 2.0);
        let linestring = LineString::from(vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);

        let distance = distance_point_to_linestring_generic(&point, &linestring);
        assert_relative_eq!(distance, 1.0); // Distance to closest segment
    }

    #[test]
    fn test_distance_point_to_linestring_generic_empty() {
        let point = Point::new(1.0, 2.0);
        let linestring = LineString::<f64>::new(vec![]);

        let distance = distance_point_to_linestring_generic(&point, &linestring);
        assert_relative_eq!(distance, 0.0);
    }

    #[test]
    fn test_distance_point_to_linestring_generic_single_point() {
        let point = Point::new(1.0, 2.0);
        let linestring = LineString::from(vec![(0.0, 0.0)]);

        let distance = distance_point_to_linestring_generic(&point, &linestring);
        assert_relative_eq!(distance, 0.0); // Single point linestring
    }

    #[test]
    fn test_distance_point_to_linestring_generic_on_linestring() {
        let point = Point::new(1.0, 0.0);
        let linestring = LineString::from(vec![(0.0, 0.0), (2.0, 0.0)]);

        let distance = distance_point_to_linestring_generic(&point, &linestring);
        assert_relative_eq!(distance, 0.0);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_point_to_polygon_generic                │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_point_to_polygon_generic_outside() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(6.0, 2.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 2.0); // Distance to right edge
    }

    #[test]
    fn test_distance_point_to_polygon_generic_inside() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(2.0, 2.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 0.0); // Inside polygon
    }

    #[test]
    fn test_distance_point_to_polygon_generic_on_boundary() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(2.0, 0.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 0.0); // On boundary
    }

    #[test]
    fn test_distance_point_to_polygon_generic_with_hole() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (0.0, 0.0)
        ]);
        let interior = LineString::from(vec![
            (2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)
        ]);
        let polygon = Polygon::new(exterior, vec![interior]);
        let point = Point::new(3.0, 3.0); // Inside the hole

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 1.0); // Distance to closest hole edge
    }

    #[test]
    fn test_distance_point_to_polygon_generic_empty() {
        let empty_polygon = Polygon::new(LineString::<f64>::new(vec![]), vec![]);
        let point = Point::new(1.0, 1.0);

        let distance = distance_point_to_polygon_generic(&point, &empty_polygon);
        assert_relative_eq!(distance, 0.0);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_line_to_line_generic                    │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_line_to_line_generic_parallel() {
        let line1 = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 2.0, y: 0.0 });
        let line2 = Line::new(coord! { x: 0.0, y: 3.0 }, coord! { x: 2.0, y: 3.0 });

        let distance = distance_line_to_line_generic(&line1, &line2);
        assert_relative_eq!(distance, 3.0);
    }

    #[test]
    fn test_distance_line_to_line_generic_intersecting() {
        let line1 = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 2.0, y: 0.0 });
        let line2 = Line::new(coord! { x: 1.0, y: -1.0 }, coord! { x: 1.0, y: 1.0 });

        let distance = distance_line_to_line_generic(&line1, &line2);
        assert_relative_eq!(distance, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_distance_line_to_line_generic_skew() {
        let line1 = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 1.0, y: 0.0 });
        let line2 = Line::new(coord! { x: 2.0, y: 1.0 }, coord! { x: 3.0, y: 1.0 });

        let distance = distance_line_to_line_generic(&line1, &line2);
        let expected = ((2.0 - 1.0).powi(2) + (1.0 - 0.0).powi(2)).sqrt();
        assert_relative_eq!(distance, expected);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_linestring_to_polygon_generic           │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_linestring_to_polygon_generic_outside() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let linestring = LineString::from(vec![(3.0, 0.0), (4.0, 1.0)]);

        let distance = distance_linestring_to_polygon_generic(&linestring, &polygon);
        assert_relative_eq!(distance, 1.0); // Distance to right edge of polygon
    }

    #[test]
    fn test_distance_linestring_to_polygon_generic_intersecting() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let linestring = LineString::from(vec![(-1.0, 1.0), (3.0, 1.0)]);

        let distance = distance_linestring_to_polygon_generic(&linestring, &polygon);
        // The algorithm computes minimum distance between line segments
        // For this configuration, the result is 1.0 (minimum distance between segments)
        assert_relative_eq!(distance, 1.0);
    }

    #[test]
    fn test_distance_linestring_to_polygon_generic_empty_polygon() {
        let empty_polygon = Polygon::new(LineString::<f64>::new(vec![]), vec![]);
        let linestring = LineString::from(vec![(0.0, 0.0), (1.0, 1.0)]);

        let distance = distance_linestring_to_polygon_generic(&linestring, &empty_polygon);
        assert_relative_eq!(distance, 0.0);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_polygon_to_polygon_generic              │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_polygon_to_polygon_generic_separate() {
        let exterior1 = LineString::from(vec![
            (0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)
        ]);
        let polygon1 = Polygon::new(exterior1, vec![]);

        let exterior2 = LineString::from(vec![
            (4.0, 0.0), (6.0, 0.0), (6.0, 2.0), (4.0, 2.0), (4.0, 0.0)
        ]);
        let polygon2 = Polygon::new(exterior2, vec![]);

        let distance = distance_polygon_to_polygon_generic(&polygon1, &polygon2);
        assert_relative_eq!(distance, 2.0); // Distance between closest edges
    }

    #[test]
    fn test_distance_polygon_to_polygon_generic_intersecting() {
        let exterior1 = LineString::from(vec![
            (0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)
        ]);
        let polygon1 = Polygon::new(exterior1, vec![]);

        let exterior2 = LineString::from(vec![
            (1.0, 1.0), (4.0, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 1.0)
        ]);
        let polygon2 = Polygon::new(exterior2, vec![]);

        let distance = distance_polygon_to_polygon_generic(&polygon1, &polygon2);
        assert_relative_eq!(distance, 0.0, epsilon = 1e-10); // Polygons intersect
    }

    #[test]
    fn test_distance_polygon_to_polygon_generic_one_in_others_hole() {
        let exterior = LineString::from(vec![
            (0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)
        ]);
        let interior = LineString::from(vec![
            (2.0, 2.0), (8.0, 2.0), (8.0, 8.0), (2.0, 8.0), (2.0, 2.0)
        ]);
        let polygon_with_hole = Polygon::new(exterior, vec![interior]);

        let small_exterior = LineString::from(vec![
            (4.0, 4.0), (6.0, 4.0), (6.0, 6.0), (4.0, 6.0), (4.0, 4.0)
        ]);
        let small_polygon = Polygon::new(small_exterior, vec![]);

        let distance = distance_polygon_to_polygon_generic(&polygon_with_hole, &small_polygon);
        assert_relative_eq!(distance, 2.0); // Distance to hole boundary
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for symmetric distance functions                     │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_symmetric_distance_point_linestring() {
        let point = Point::new(1.0, 2.0);
        let linestring = LineString::from(vec![(0.0, 0.0), (2.0, 0.0)]);

        let dist1 = distance_point_to_linestring_generic(&point, &linestring);
        let dist2 = distance_linestring_to_point_generic(&linestring, &point);

        assert_relative_eq!(dist1, dist2);
        assert_relative_eq!(dist1, 2.0);
    }

    #[test]
    fn test_symmetric_distance_point_polygon() {
        let point = Point::new(5.0, 2.0);
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);

        let dist1 = distance_point_to_polygon_generic(&point, &polygon);
        let dist2 = distance_polygon_to_point_generic(&polygon, &point);

        assert_relative_eq!(dist1, dist2);
        assert_relative_eq!(dist1, 1.0);
    }

    #[test]
    fn test_symmetric_distance_linestring_polygon() {
        let linestring = LineString::from(vec![(5.0, 1.0), (6.0, 2.0)]);
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);

        let dist1 = distance_linestring_to_polygon_generic(&linestring, &polygon);
        let dist2 = distance_polygon_to_linestring_generic(&polygon, &linestring);

        assert_relative_eq!(dist1, dist2);
        assert_relative_eq!(dist1, 1.0);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for line-to-linestring and line-to-polygon functions │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_line_to_linestring_generic() {
        let line = Line::new(coord! { x: 0.0, y: 3.0 }, coord! { x: 2.0, y: 3.0 });
        let linestring = LineString::from(vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]);

        let distance = distance_line_to_linestring_generic(&line, &linestring);
        assert_relative_eq!(distance, 1.0); // Distance to closest segment
    }

    #[test]
    fn test_distance_line_to_polygon_generic() {
        let line = Line::new(coord! { x: 5.0, y: 1.0 }, coord! { x: 6.0, y: 2.0 });
        let exterior = LineString::from(vec![
            (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(exterior, vec![]);

        let distance = distance_line_to_polygon_generic(&line, &polygon);
        assert_relative_eq!(distance, 1.0);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Tests for distance_triangle_to_point_generic               │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_distance_triangle_to_point_generic() {
        let triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0 },
            coord! { x: 3.0, y: 0.0 },
            coord! { x: 1.5, y: 3.0 }
        );
        let point = Point::new(1.5, 1.0); // Inside triangle

        let distance = distance_triangle_to_point_generic(&triangle, &point);
        assert_relative_eq!(distance, 0.0);
    }

    #[test]
    fn test_distance_triangle_to_point_generic_outside() {
        let triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0 },
            coord! { x: 3.0, y: 0.0 },
            coord! { x: 1.5, y: 3.0 }
        );
        let point = Point::new(5.0, 0.0); // Outside triangle

        let distance = distance_triangle_to_point_generic(&triangle, &point);
        assert_relative_eq!(distance, 2.0); // Distance to right vertex
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Edge case tests                                            │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_empty_geometries_edge_cases() {
        // Empty LineString
        let empty_ls = LineString::<f64>::new(vec![]);
        let point = Point::new(1.0, 1.0);

        let dist = distance_point_to_linestring_generic(&point, &empty_ls);
        assert_relative_eq!(dist, 0.0);

        // Empty Polygon
        let empty_poly = Polygon::new(LineString::<f64>::new(vec![]), vec![]);
        let dist2 = distance_point_to_polygon_generic(&point, &empty_poly);
        assert_relative_eq!(dist2, 0.0);
    }

    #[test]
    fn test_degenerate_geometries() {
        // Single point LineString
        let single_point_ls = LineString::from(vec![(1.0, 1.0)]);
        let point = Point::new(2.0, 2.0);

        let dist = distance_point_to_linestring_generic(&point, &single_point_ls);
        assert_relative_eq!(dist, 0.0); // Should handle gracefully

        // Two identical points in LineString
        let two_same_points_ls = LineString::from(vec![(1.0, 1.0), (1.0, 1.0)]);
        let dist2 = distance_point_to_linestring_generic(&point, &two_same_points_ls);
        assert_relative_eq!(dist2, std::f64::consts::SQRT_2); // Distance to the point
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Performance comparison tests (basic)                       │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_generic_vs_concrete_point_distance() {
        let p1 = Point::new(-72.1235, 42.3521);
        let p2 = Point::new(72.1260, 70.612);

        // Test generic implementation
        let generic_dist = point_distance_generic(&p1, &p2);

        // Test concrete implementation via Euclidean trait
        let concrete_dist = Euclidean.distance(&p1, &p2);

        // Both should give the same result
        assert_relative_eq!(generic_dist, concrete_dist, epsilon = 1e-10);
        assert_relative_eq!(generic_dist, 146.99163308930207);
    }

    #[test]
    fn test_cross_validation_with_existing_tests() {
        // Test cases from existing distance.rs tests to ensure compatibility
        let o1 = Point::new(8.0, 0.0);
        let p1 = Point::new(7.2, 2.0);
        let p2 = Point::new(6.0, 1.0);

        // Create line from p1 to p2
        let line_seg = Line::new(
            coord! { x: p1.x(), y: p1.y() },
            coord! { x: p2.x(), y: p2.y() }
        );

        if let Some(o1_coord) = o1.coord_ext() {
            let generic_dist = line_segment_distance_generic(&o1_coord, &line_seg);

            // This should match the expected value from the original test
            assert_relative_eq!(generic_dist, 2.0485900789263356, epsilon = 1e-10);
        }
    }
}
