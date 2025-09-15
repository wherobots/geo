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

pub fn nearest_neighbour_distance<F: GeoFloat>(geom1: &LineString<F>, geom2: &LineString<F>) -> F {
    let tree_a = RTree::bulk_load(geom1.lines().map(CachedEnvelope::new).collect());
    let tree_b = RTree::bulk_load(geom2.lines().map(CachedEnvelope::new).collect());

    let mut min_distance: F = Bounded::max_value();

    for line1 in geom1.lines() {
        for line2 in geom2.lines() {
            let line_distance = distance_line_to_line_generic(&line1, &line2);
            min_distance = min_distance.min(line_distance);

            // Early exit if we found an intersection
            if line_distance == F::zero() {
                return F::zero();
            }
        }
    }

    let point_line_dist = geom2
        .points()
        .fold(Bounded::max_value(), |acc: F, point| {
            let nearest = tree_a.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        })
        .min(geom1.points().fold(Bounded::max_value(), |acc, point| {
            let nearest = tree_b.nearest_neighbor(&point).unwrap();
            acc.min(Euclidean.distance(nearest as &Line<F>, &point))
        }));

    min_distance.min(point_line_dist)
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

        // Calculate distance to exterior ring using proper line-to-line distance
        for line1 in linestring.lines() {
            for line2 in exterior.lines() {
                let line_dist = distance_line_to_line_generic(&line1, &line2);
                min_dist = min_dist.min(line_dist);

                // Early exit if we found an intersection
                if line_dist == F::zero() {
                    return F::zero();
                }
            }
        }

        // Also calculate distance to interior rings (holes)
        for interior in polygon.interiors_ext() {
            for line1 in linestring.lines() {
                for line2 in interior.lines() {
                    let line_dist = distance_line_to_line_generic(&line1, &line2);
                    min_dist = min_dist.min(line_dist);

                    // Early exit if we found an intersection
                    if line_dist == F::zero() {
                        return F::zero();
                    }
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
        // The linestrings intersect at (2,0), so distance should be 0.0
        assert_relative_eq!(distance, 0.0);
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
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let coord = coord! { x: 2.0, y: 2.0 };

        assert!(ring_contains_coord(&ring, coord));
    }

    #[test]
    fn test_ring_contains_coord_outside() {
        let ring = LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let coord = coord! { x: 5.0, y: 2.0 };

        assert!(!ring_contains_coord(&ring, coord));
    }

    #[test]
    fn test_ring_contains_coord_on_boundary() {
        let ring = LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let coord = coord! { x: 2.0, y: 0.0 };

        assert!(!ring_contains_coord(&ring, coord)); // On boundary = false
    }

    #[test]
    fn test_ring_contains_coord_triangle() {
        let ring = LineString::from(vec![(0.0, 0.0), (3.0, 0.0), (1.5, 2.0), (0.0, 0.0)]);
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
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(6.0, 2.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 2.0); // Distance to right edge
    }

    #[test]
    fn test_distance_point_to_polygon_generic_inside() {
        let exterior = LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(2.0, 2.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 0.0); // Inside polygon
    }

    #[test]
    fn test_distance_point_to_polygon_generic_on_boundary() {
        let exterior = LineString::from(vec![
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let point = Point::new(2.0, 0.0);

        let distance = distance_point_to_polygon_generic(&point, &polygon);
        assert_relative_eq!(distance, 0.0); // On boundary
    }

    #[test]
    fn test_distance_point_to_polygon_generic_with_hole() {
        let exterior = LineString::from(vec![
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 6.0),
            (0.0, 6.0),
            (0.0, 0.0),
        ]);
        let interior = LineString::from(vec![
            (2.0, 2.0),
            (4.0, 2.0),
            (4.0, 4.0),
            (2.0, 4.0),
            (2.0, 2.0),
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
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let linestring = LineString::from(vec![(3.0, 0.0), (4.0, 1.0)]);

        let distance = distance_linestring_to_polygon_generic(&linestring, &polygon);
        assert_relative_eq!(distance, 1.0); // Distance to right edge of polygon
    }

    #[test]
    fn test_distance_linestring_to_polygon_generic_intersecting() {
        let exterior = LineString::from(vec![
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]);
        let polygon = Polygon::new(exterior, vec![]);
        let linestring = LineString::from(vec![(-1.0, 1.0), (3.0, 1.0)]);

        let distance = distance_linestring_to_polygon_generic(&linestring, &polygon);
        // The linestring intersects the polygon, so distance should be 0.0
        assert_relative_eq!(distance, 0.0);
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
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (0.0, 0.0),
        ]);
        let polygon1 = Polygon::new(exterior1, vec![]);

        let exterior2 = LineString::from(vec![
            (4.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (4.0, 2.0),
            (4.0, 0.0),
        ]);
        let polygon2 = Polygon::new(exterior2, vec![]);

        let distance = distance_polygon_to_polygon_generic(&polygon1, &polygon2);
        assert_relative_eq!(distance, 2.0); // Distance between closest edges
    }

    #[test]
    fn test_distance_polygon_to_polygon_generic_intersecting() {
        let exterior1 = LineString::from(vec![
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (0.0, 3.0),
            (0.0, 0.0),
        ]);
        let polygon1 = Polygon::new(exterior1, vec![]);

        let exterior2 = LineString::from(vec![
            (1.0, 1.0),
            (4.0, 1.0),
            (4.0, 4.0),
            (1.0, 4.0),
            (1.0, 1.0),
        ]);
        let polygon2 = Polygon::new(exterior2, vec![]);

        let distance = distance_polygon_to_polygon_generic(&polygon1, &polygon2);
        assert_relative_eq!(distance, 0.0, epsilon = 1e-10); // Polygons intersect
    }

    #[test]
    fn test_distance_polygon_to_polygon_generic_one_in_others_hole() {
        let exterior = LineString::from(vec![
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]);
        let interior = LineString::from(vec![
            (2.0, 2.0),
            (8.0, 2.0),
            (8.0, 8.0),
            (2.0, 8.0),
            (2.0, 2.0),
        ]);
        let polygon_with_hole = Polygon::new(exterior, vec![interior]);

        let small_exterior = LineString::from(vec![
            (4.0, 4.0),
            (6.0, 4.0),
            (6.0, 6.0),
            (4.0, 6.0),
            (4.0, 4.0),
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
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
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
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
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
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (0.0, 0.0),
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
            coord! { x: 1.5, y: 3.0 },
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
            coord! { x: 1.5, y: 3.0 },
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
            coord! { x: p2.x(), y: p2.y() },
        );

        if let Some(o1_coord) = o1.coord_ext() {
            let generic_dist = line_segment_distance_generic(&o1_coord, &line_seg);

            // This should match the expected value from the original test
            assert_relative_eq!(generic_dist, 2.0485900789263356, epsilon = 1e-10);
        }
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Property-based tests with random inputs                    │
    // └────────────────────────────────────────────────────────────┘

    fn generate_random_point(seed: u64) -> Point<f64> {
        // Simple LCG for deterministic "random" numbers
        let mut rng = seed;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let x = ((rng >> 16) as i16) as f64 * 0.001;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let y = ((rng >> 16) as i16) as f64 * 0.001;
        Point::new(x, y)
    }

    fn generate_random_line(seed: u64) -> Line<f64> {
        let mut rng = seed;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let x1 = ((rng >> 16) as i16) as f64 * 0.001;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let y1 = ((rng >> 16) as i16) as f64 * 0.001;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let x2 = ((rng >> 16) as i16) as f64 * 0.001;
        rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
        let y2 = ((rng >> 16) as i16) as f64 * 0.001;
        Line::new(coord! { x: x1, y: y1 }, coord! { x: x2, y: y2 })
    }

    fn generate_random_linestring(seed: u64, num_points: usize) -> LineString<f64> {
        let mut rng = seed;
        let mut points = Vec::new();
        for _ in 0..num_points {
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            let x = ((rng >> 16) as i16) as f64 * 0.001;
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            let y = ((rng >> 16) as i16) as f64 * 0.001;
            points.push((x, y));
        }
        LineString::from(points)
    }

    fn generate_random_polygon(seed: u64, num_exterior_points: usize) -> Polygon<f64> {
        let mut rng = seed;
        let mut points = Vec::new();

        // Generate points around a circle to ensure a valid polygon
        let center_x = 0.0;
        let center_y = 0.0;
        let radius = 10.0;

        for i in 0..num_exterior_points {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / num_exterior_points as f64;
            // Add some random noise
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            let noise = ((rng >> 16) as i16) as f64 * 0.0001;
            let x = center_x + (radius + noise) * angle.cos();
            let y = center_y + (radius + noise) * angle.sin();
            points.push((x, y));
        }

        // Close the polygon
        if !points.is_empty() {
            points.push(points[0]);
        }

        Polygon::new(LineString::from(points), vec![])
    }

    #[test]
    fn test_random_point_to_point_distance() {
        // Test point-to-point distance with random inputs
        for i in 0..100 {
            let seed1 = 12345 + i * 17;
            let seed2 = 54321 + i * 23;

            let p1 = generate_random_point(seed1);
            let p2 = generate_random_point(seed2);

            let concrete_dist = Euclidean.distance(&p1, &p2);
            let generic_dist = point_distance_generic(&p1, &p2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-12,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn test_random_point_to_linestring_distance() {
        // Test point-to-linestring distance with random inputs
        for i in 0..100 {
            let seed1 = 11111 + i * 31;
            let seed2 = 22222 + i * 37;

            let point = generate_random_point(seed1);
            let linestring = generate_random_linestring(seed2, 3 + (i % 5) as usize); // 3-7 points

            let concrete_dist = Euclidean.distance(&point, &linestring);
            let generic_dist = distance_point_to_linestring_generic(&point, &linestring);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-12,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn test_random_point_to_polygon_distance() {
        // Test point-to-polygon distance with random inputs
        for i in 0..100 {
            let seed1 = 33333 + i * 41;
            let seed2 = 44444 + i * 43;

            let point = generate_random_point(seed1);
            let polygon = generate_random_polygon(seed2, 4 + (i % 4) as usize); // 4-7 sides

            let concrete_dist = Euclidean.distance(&point, &polygon);
            let generic_dist = distance_point_to_polygon_generic(&point, &polygon);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-10,
                max_relative = 1e-10
            );
        }
    }

    #[test]
    fn test_random_line_to_line_distance() {
        // Test line-to-line distance with random inputs
        for i in 0..100 {
            let seed1 = 55555 + i * 47;
            let seed2 = 66666 + i * 53;

            let line1 = generate_random_line(seed1);
            let line2 = generate_random_line(seed2);

            let concrete_dist = Euclidean.distance(&line1, &line2);
            let generic_dist = distance_line_to_line_generic(&line1, &line2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-12,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn test_random_linestring_to_linestring_distance() {
        // Test linestring-to-linestring distance with random inputs
        for i in 0..100 {
            let seed1 = 77777 + i * 59;
            let seed2 = 88888 + i * 61;

            let ls1 = generate_random_linestring(seed1, 3 + (i % 3) as usize); // 3-5 points
            let ls2 = generate_random_linestring(seed2, 3 + ((i + 1) % 3) as usize); // 3-5 points

            let concrete_dist = Euclidean.distance(&ls1, &ls2);
            // Use our actual generic implementation via nearest_neighbour_distance
            let generic_dist = nearest_neighbour_distance(&ls1, &ls2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-10,
                max_relative = 1e-10
            );
        }
    }

    #[test]
    fn test_random_polygon_to_polygon_distance() {
        // Test polygon-to-polygon distance with random inputs
        for i in 0..100 {
            let seed1 = 99999 + i * 67;
            let seed2 = 10101 + i * 71;

            let poly1 = generate_random_polygon(seed1, 4 + (i % 3) as usize); // 4-6 sides
            let poly2 = generate_random_polygon(seed2, 4 + ((i + 1) % 3) as usize); // 4-6 sides

            let concrete_dist = Euclidean.distance(&poly1, &poly2);
            let generic_dist = distance_polygon_to_polygon_generic(&poly1, &poly2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-8,
                max_relative = 1e-8
            );
        }
    }

    #[test]
    fn test_random_line_to_polygon_distance() {
        // Test line-to-polygon distance with random inputs
        for i in 0..100 {
            let seed1 = 12121 + i * 73;
            let seed2 = 13131 + i * 79;

            let line = generate_random_line(seed1);
            let polygon = generate_random_polygon(seed2, 4 + (i % 3) as usize); // 4-6 sides

            let concrete_dist = Euclidean.distance(&line, &polygon);
            let generic_dist = distance_line_to_polygon_generic(&line, &polygon);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-10,
                max_relative = 1e-10
            );
        }
    }

    #[test]
    fn test_random_linestring_to_polygon_distance() {
        // Test linestring-to-polygon distance with random inputs
        for i in 0..100 {
            let seed1 = 14141 + i * 83;
            let seed2 = 15151 + i * 89;

            let linestring = generate_random_linestring(seed1, 3 + (i % 3) as usize); // 3-5 points
            let polygon = generate_random_polygon(seed2, 4 + (i % 3) as usize); // 4-6 sides

            let concrete_dist = Euclidean.distance(&linestring, &polygon);
            let generic_dist = distance_linestring_to_polygon_generic(&linestring, &polygon);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-8,
                max_relative = 1e-8
            );
        }
    }

    #[test]
    fn test_random_symmetry_properties() {
        // Test symmetry properties with random inputs
        for i in 0..100 {
            let seed1 = 16161 + i * 97;
            let seed2 = 17171 + i * 101;

            let point = generate_random_point(seed1);
            let linestring = generate_random_linestring(seed2, 4);

            // Test point-linestring symmetry
            let dist1 = distance_point_to_linestring_generic(&point, &linestring);
            let dist2 = distance_linestring_to_point_generic(&linestring, &point);
            assert_relative_eq!(dist1, dist2, epsilon = 1e-12);

            // Test with polygon
            if i % 2 == 0 {
                let polygon = generate_random_polygon(seed1 + seed2, 5);
                let dist3 = distance_point_to_polygon_generic(&point, &polygon);
                let dist4 = distance_polygon_to_point_generic(&polygon, &point);
                assert_relative_eq!(dist3, dist4, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_random_edge_cases_and_boundaries() {
        // Test edge cases with specific patterns
        for i in 0..100 {
            // Same point distance should be zero
            let point = generate_random_point(12345 + i);
            let same_point_dist = point_distance_generic(&point, &point);
            assert_relative_eq!(same_point_dist, 0.0);

            // Zero-length line segment
            let coord = coord! { x: point.x(), y: point.y() };
            let zero_line = Line::new(coord, coord);
            let dist_to_zero_line = line_segment_distance_generic(&coord, &zero_line);
            assert_relative_eq!(dist_to_zero_line, 0.0);

            // Point on line segment should have zero distance
            let seed = 54321 + i * 13;
            let line = generate_random_line(seed);
            let start_coord = line.start_coord();
            let dist_to_start = line_segment_distance_generic(&start_coord, &line);
            assert_relative_eq!(dist_to_start, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_random_large_coordinates() {
        // Test with large coordinate values to check numerical stability
        for i in 0..100 {
            let mut rng: u64 = 98765 + i * 107;

            // Generate large coordinates
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            let scale = 1e6 + (rng % 1000000) as f64;

            let p1 = Point::new(scale, scale * 0.5);
            let p2 = Point::new(scale * 1.1, scale * 0.7);

            let concrete_dist = Euclidean.distance(&p1, &p2);
            let generic_dist = point_distance_generic(&p1, &p2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-10,
                max_relative = 1e-10
            );
        }
    }

    #[test]
    fn test_random_small_coordinates() {
        // Test with very small coordinate values to check numerical precision
        for i in 0..100 {
            let mut rng: u64 = 13579 + i * 109;

            // Generate small coordinates
            rng = rng.wrapping_mul(1103515245).wrapping_add(12345);
            let scale = 1e-6 * (1.0 + (rng % 100) as f64 * 0.01);

            let p1 = Point::new(scale, scale * 0.5);
            let p2 = Point::new(scale * 1.1, scale * 0.7);

            let concrete_dist = Euclidean.distance(&p1, &p2);
            let generic_dist = point_distance_generic(&p1, &p2);

            assert_relative_eq!(
                concrete_dist,
                generic_dist,
                epsilon = 1e-15,
                max_relative = 1e-12
            );
        }
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Geometric Edge Cases Tests                                 │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_collinear_linestring_geometries() {
        // Test linestrings where all points are collinear
        let collinear_ls1 = LineString::from(vec![
            (0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)
        ]);
        let collinear_ls2 = LineString::from(vec![
            (0.0, 1.0), (1.0, 2.0), (2.0, 3.0)
        ]);

        let concrete_dist = Euclidean.distance(&collinear_ls1, &collinear_ls2);
        let generic_dist = nearest_neighbour_distance(&collinear_ls1, &collinear_ls2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-10);
        // Distance should be sqrt(2)/2 (perpendicular distance between parallel lines)
        assert_relative_eq!(concrete_dist, std::f64::consts::SQRT_2 / 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_degenerate_triangle_as_line() {
        // Triangle where all three points are collinear (degenerate triangle)
        let degenerate_triangle = Triangle::new(
            coord! { x: 0.0, y: 0.0 },
            coord! { x: 1.0, y: 1.0 },
            coord! { x: 2.0, y: 2.0 },
        );
        let point = Point::new(0.0, 1.0);

        let concrete_dist = Euclidean.distance(&degenerate_triangle, &point);
        let generic_dist = distance_triangle_to_point_generic(&degenerate_triangle, &point);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-10);
        // Distance should be sqrt(2)/2 (distance from point to line y=x)
        assert_relative_eq!(concrete_dist, std::f64::consts::SQRT_2 / 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_self_intersecting_polygon() {
        // Create a bowtie/figure-8 shaped self-intersecting polygon
        let self_intersecting = LineString::from(vec![
            (0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0), (0.0, 0.0)
        ]);
        let polygon = Polygon::new(self_intersecting, vec![]);
        let point = Point::new(3.0, 1.0); // Outside the polygon

        let concrete_dist = Euclidean.distance(&point, &polygon);
        let generic_dist = distance_point_to_polygon_generic(&point, &polygon);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-10);
        assert_relative_eq!(concrete_dist, 1.0, epsilon = 1e-10); // Distance to closest edge
    }

    #[test]
    fn test_nearly_touching_geometries() {
        // Test geometries separated by very small distances
        let epsilon_dist = 1e-12;

        let line1 = Line::new(coord! { x: 0.0, y: 0.0 }, coord! { x: 1.0, y: 0.0 });
        let line2 = Line::new(coord! { x: 0.0, y: epsilon_dist }, coord! { x: 1.0, y: epsilon_dist });

        let concrete_dist = Euclidean.distance(&line1, &line2);
        let generic_dist = distance_line_to_line_generic(&line1, &line2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-15);
        assert_relative_eq!(concrete_dist, epsilon_dist, epsilon = 1e-15);
    }

    #[test]
    fn test_very_close_but_separate_polygons() {
        // Two polygons separated by extremely small distance
        let tiny_gap = 1e-14;

        let poly1_exterior = LineString::from(vec![
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)
        ]);
        let poly1 = Polygon::new(poly1_exterior, vec![]);

        let poly2_exterior = LineString::from(vec![
            (1.0 + tiny_gap, 0.0), (2.0 + tiny_gap, 0.0),
            (2.0 + tiny_gap, 1.0), (1.0 + tiny_gap, 1.0), (1.0 + tiny_gap, 0.0)
        ]);
        let poly2 = Polygon::new(poly2_exterior, vec![]);

        let concrete_dist = Euclidean.distance(&poly1, &poly2);
        let generic_dist = distance_polygon_to_polygon_generic(&poly1, &poly2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-15);
        assert_relative_eq!(concrete_dist, tiny_gap, epsilon = 1e-16);
    }

    #[test]
    fn test_overlapping_but_not_intersecting_linestrings() {
        // LineStrings that overlap in projection but are at different heights
        let ls1 = LineString::from(vec![(0.0, 0.0), (2.0, 0.0)]);
        let ls2 = LineString::from(vec![(1.0, 1e-13), (3.0, 1e-13)]);

        let concrete_dist = Euclidean.distance(&ls1, &ls2);
        let generic_dist = nearest_neighbour_distance(&ls1, &ls2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-15);
        assert_relative_eq!(concrete_dist, 1e-13, epsilon = 1e-16);
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Numerical Precision Tests                                  │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_very_close_but_non_zero_distances() {
        // Test extremely small but non-zero distances to check floating-point precision
        let test_cases = [1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10];

        for &tiny_dist in &test_cases {
            let p1 = Point::new(0.0, 0.0);
            let p2 = Point::new(tiny_dist, 0.0);

            let concrete_dist = Euclidean.distance(&p1, &p2);
            let generic_dist = point_distance_generic(&p1, &p2);

            assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-16);
            assert_relative_eq!(concrete_dist, tiny_dist, epsilon = 1e-16);
            assert!(concrete_dist > 0.0, "Distance should be positive for tiny_dist = {}", tiny_dist);
        }
    }

    #[test]
    fn test_numerical_precision_near_floating_point_limits() {
        // Test with coordinates that produce distances near floating-point precision limits
        let base = 1.0;
        let tiny_offset = std::f64::EPSILON * 10.0; // Slightly above machine epsilon

        let p1 = Point::new(base, base);
        let p2 = Point::new(base + tiny_offset, base);

        let concrete_dist = Euclidean.distance(&p1, &p2);
        let generic_dist = point_distance_generic(&p1, &p2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-15);
        assert!(concrete_dist > 0.0);
        assert!(concrete_dist < 1e-14); // Should be very small but measurable
    }

    #[test]
    fn test_precision_with_large_coordinate_differences() {
        // Test with one geometry having small coordinates and another having large coordinates
        let small_point = Point::new(1e-10, 1e-10);
        let large_polygon = Polygon::new(
            LineString::from(vec![
                (1e8, 1e8), (1e8 + 1.0, 1e8), (1e8 + 1.0, 1e8 + 1.0), (1e8, 1e8 + 1.0), (1e8, 1e8)
            ]),
            vec![]
        );

        let concrete_dist = Euclidean.distance(&small_point, &large_polygon);
        let generic_dist = distance_point_to_polygon_generic(&small_point, &large_polygon);

        assert_relative_eq!(concrete_dist, generic_dist, max_relative = 1e-10);
        assert!(concrete_dist > 1e7); // Should be very large distance
    }

    // ┌────────────────────────────────────────────────────────────┐
    // │ Robustness Tests                                           │
    // └────────────────────────────────────────────────────────────┘

    #[test]
    fn test_nan_coordinate_handling() {
        // Test behavior with NaN coordinates
        let nan_point = Point::new(f64::NAN, 0.0);
        let normal_point = Point::new(1.0, 1.0);

        let distance = point_distance_generic(&nan_point, &normal_point);

        // Distance involving NaN should be NaN
        assert!(distance.is_nan(), "Distance with NaN coordinate should be NaN");
    }

    #[test]
    fn test_infinity_coordinate_handling() {
        // Test behavior with infinite coordinates
        let inf_point = Point::new(f64::INFINITY, 0.0);
        let normal_point = Point::new(1.0, 1.0);

        let distance = point_distance_generic(&inf_point, &normal_point);

        // Distance involving infinity should be infinity
        assert!(distance.is_infinite(), "Distance with infinite coordinate should be infinite");
    }

    #[test]
    fn test_negative_infinity_coordinate_handling() {
        // Test behavior with negative infinite coordinates
        let neg_inf_point = Point::new(f64::NEG_INFINITY, 0.0);
        let normal_point = Point::new(1.0, 1.0);

        let distance = point_distance_generic(&neg_inf_point, &normal_point);

        // Distance involving negative infinity should be infinity
        assert!(distance.is_infinite(), "Distance with negative infinite coordinate should be infinite");
    }

    #[test]
    fn test_mixed_special_values() {
        // Test combinations of NaN and infinity
        let nan_point = Point::new(f64::NAN, f64::INFINITY);
        let inf_point = Point::new(f64::INFINITY, f64::NEG_INFINITY);

        let distance = point_distance_generic(&nan_point, &inf_point);

        // Any operation involving NaN should result in NaN or Infinity depending on the math
        // Since we're using hypot which can handle NaN differently, let's test that it's either NaN or infinite
        assert!(distance.is_nan() || distance.is_infinite(),
                "Distance involving NaN and Infinity should be NaN or Infinite, got: {}", distance);
    }

    #[test]
    fn test_subnormal_number_handling() {
        // Test with subnormal (denormalized) numbers
        let subnormal = f64::MIN_POSITIVE / 2.0; // This creates a subnormal number
        assert!(subnormal > 0.0 && subnormal < f64::MIN_POSITIVE);

        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(subnormal, 0.0);

        let concrete_dist = Euclidean.distance(&p1, &p2);
        let generic_dist = point_distance_generic(&p1, &p2);

        assert_relative_eq!(concrete_dist, generic_dist, epsilon = 1e-16);
        assert_relative_eq!(concrete_dist, subnormal, epsilon = 1e-16);
        assert!(concrete_dist > 0.0);
    }

    #[test]
    fn test_zero_vs_negative_zero() {
        // Test behavior with positive zero vs negative zero
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(-0.0, -0.0); // Negative zero

        let distance = point_distance_generic(&p1, &p2);

        // Distance between +0 and -0 should be exactly 0
        assert_eq!(distance, 0.0, "Distance between +0 and -0 should be exactly 0");
    }
}
