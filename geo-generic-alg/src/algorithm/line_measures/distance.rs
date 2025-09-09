use crate::{CoordFloat, Point};
use geo_traits::{CoordTrait, LineStringTrait, PointTrait};
use geo_traits_ext::*;

/// Calculate the minimum distance between two geometries.
pub trait Distance<F, Origin, Destination> {
    /// Note that not all implementations support all geometry combinations, but at least `Point` to `Point`
    /// is supported.
    /// See [specific implementations](#implementors) for details.
    ///
    /// # Units
    ///
    /// - `origin`, `destination`: geometry where the units of x/y depend on the trait implementation.
    /// - returns: depends on the trait implementation.
    ///
    /// # Examples
    ///
    /// ```
    /// use geo::{Haversine, Euclidean, Distance, Point};
    /// let p1: Point = Point::new(0.0, 0.0);
    /// let p2: Point = Point::new(0.0, 2.0);
    ///
    /// assert_eq!(Euclidean.distance(p1, p2), 2.0);
    ///
    /// // The units of the output depend on the metric space.
    /// // In the case of [`Haversine`], it's meters.
    /// // See the documentation for each metric space for details.
    /// assert_eq!(Haversine.distance(p1, p2).round(), 222_390.0);
    /// ```
    fn distance(&self, origin: Origin, destination: Destination) -> F;
}

/// Extension trait that enables generic distance calculations for WKB and other generic geometry types.
///
/// This provides WKB-compatible distance calculations that work with any geometry type
/// implementing the geo-traits-ext pattern, enabling distance calculations for
/// WKB geometries, geoarrow, and other generic geometry representations.
///
/// # Examples
/// ```
/// use geo_generic_alg::algorithm::line_measures::{DistanceExt, Euclidean};
/// // Two WKB Points: (0,0) and (3,4)
/// let wkb_bytes1: &[u8] = &[
///     1, // little endian
///     1, 0, 0, 0, // geometry type: Point (1)
///     0, 0, 0, 0, 0, 0, 0, 0, // x = 0.0
///     0, 0, 0, 0, 0, 0, 0, 0, // y = 0.0
/// ];
/// let wkb_bytes2: &[u8] = &[
///     1, // little endian
///     1, 0, 0, 0, // geometry type: Point (1)
///     0, 0, 0, 0, 0, 0, 8, 64, // x = 3.0
///     0, 0, 0, 0, 0, 0, 16, 64, // y = 4.0
/// ];
/// let wkb_geom1 = geo_generic_tests::wkb::reader::read_wkb(wkb_bytes1).unwrap();
/// let wkb_geom2 = geo_generic_tests::wkb::reader::read_wkb(wkb_bytes2).unwrap();
///
/// let distance = wkb_geom1.distance_ext(&Euclidean, &wkb_geom2);
/// assert_eq!(distance, 5.0);
/// ```
pub trait DistanceExt<F: CoordFloat> {
    /// Calculate the distance to another geometry using the given metric space.
    fn distance_ext(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>, other: &Self) -> F;
}

// Implementation for WKB and other generic geometries using the type-tag pattern
impl<F, G> DistanceExt<F> for G
where
    F: CoordFloat,
    G: GeoTraitExtWithTypeTag + DistanceTrait<F, G::Tag>,
{
    fn distance_ext(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>, other: &Self) -> F {
        self.distance_trait(metric_space, other)
    }
}

// Internal trait that handles the actual distance computation for different geometry types
trait DistanceTrait<F, GT: GeoTypeTag>
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F;
}

// Implementation for Point geometries
impl<F, P: PointTraitExt<T = F>> DistanceTrait<F, PointTag> for P
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        if let (Some(coord1), Some(coord2)) = (self.coord(), other.coord()) {
            let p1 = Point::new(coord1.x(), coord1.y());
            let p2 = Point::new(coord2.x(), coord2.y());
            metric_space.distance(p1, p2)
        } else {
            F::zero()
        }
    }
}

// Implementation for Line geometries
impl<F, L: LineTraitExt<T = F>> DistanceTrait<F, LineTag> for L
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // For Line-to-Line distance, we find the minimum distance between the four endpoint pairs
        let line1_start = Point::new(self.start_coord().x, self.start_coord().y);
        let line1_end = Point::new(self.end_coord().x, self.end_coord().y);
        let line2_start = Point::new(other.start_coord().x, other.start_coord().y);
        let line2_end = Point::new(other.end_coord().x, other.end_coord().y);

        // Calculate distances between all endpoint combinations
        let d1 = metric_space.distance(line1_start, line2_start);
        let d2 = metric_space.distance(line1_start, line2_end);
        let d3 = metric_space.distance(line1_end, line2_start);
        let d4 = metric_space.distance(line1_end, line2_end);

        // Return minimum distance
        d1.min(d2).min(d3).min(d4)
    }
}

// Implementation for LineString geometries
impl<F, LS: LineStringTraitExt<T = F>> DistanceTrait<F, LineStringTag> for LS
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // For LineString-to-LineString distance, find minimum distance between all point pairs
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        for coord1 in self.coords_ext() {
            let p1 = Point::new(coord1.x(), coord1.y());
            for coord2 in other.coords_ext() {
                let p2 = Point::new(coord2.x(), coord2.y());
                let dist = metric_space.distance(p1, p2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// Implementation for Polygon geometries
impl<F, P: PolygonTraitExt<T = F>> DistanceTrait<F, PolygonTag> for P
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // For Polygon-to-Polygon distance, find minimum distance between all exterior ring vertices
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        if let (Some(exterior1), Some(exterior2)) = (self.exterior(), other.exterior()) {
            for coord1 in exterior1.coords() {
                let p1 = Point::new(coord1.x(), coord1.y());
                for coord2 in exterior2.coords() {
                    let p2 = Point::new(coord2.x(), coord2.y());
                    let dist = metric_space.distance(p1, p2);
                    if dist < min_distance {
                        min_distance = dist;
                    }
                }
            }
        }

        min_distance
    }
}

// For MultiPoint geometries
impl<F, MP: MultiPointTraitExt<T = F>> DistanceTrait<F, MultiPointTag> for MP
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // Find minimum distance between all point pairs from both MultiPoints
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        for point1 in self.points_ext() {
            if let Some(coord1) = point1.coord() {
                let p1 = Point::new(coord1.x(), coord1.y());
                for point2 in other.points_ext() {
                    if let Some(coord2) = point2.coord() {
                        let p2 = Point::new(coord2.x(), coord2.y());
                        let dist = metric_space.distance(p1, p2);
                        if dist < min_distance {
                            min_distance = dist;
                        }
                    }
                }
            }
        }

        min_distance
    }
}

// For MultiLineString geometries
impl<F, MLS: MultiLineStringTraitExt<T = F>> DistanceTrait<F, MultiLineStringTag> for MLS
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // Find minimum distance between all LineString pairs from both MultiLineStrings
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        for ls1 in self.line_strings_ext() {
            for ls2 in other.line_strings_ext() {
                let dist = ls1.distance_trait(metric_space, &ls2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// For MultiPolygon geometries
impl<F, MP: MultiPolygonTraitExt<T = F>> DistanceTrait<F, MultiPolygonTag> for MP
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // Find minimum distance between all Polygon pairs from both MultiPolygons
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        for poly1 in self.polygons_ext() {
            for poly2 in other.polygons_ext() {
                let dist = poly1.distance_trait(metric_space, &poly2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// For Rectangle geometries
impl<F, R: RectTraitExt<T = F>> DistanceTrait<F, RectTag> for R
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // For Rectangle-to-Rectangle distance, find minimum distance between corner points
        let corners1 = [
            Point::new(self.min().x(), self.min().y()), // bottom-left
            Point::new(self.max().x(), self.min().y()), // bottom-right
            Point::new(self.max().x(), self.max().y()), // top-right
            Point::new(self.min().x(), self.max().y()), // top-left
        ];

        let corners2 = [
            Point::new(other.min().x(), other.min().y()),
            Point::new(other.max().x(), other.min().y()),
            Point::new(other.max().x(), other.max().y()),
            Point::new(other.min().x(), other.max().y()),
        ];

        let mut min_distance = F::from(f64::INFINITY).unwrap();
        for &p1 in &corners1 {
            for &p2 in &corners2 {
                let dist = metric_space.distance(p1, p2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// For Triangle geometries
impl<F, T: TriangleTraitExt<T = F>> DistanceTrait<F, TriangleTag> for T
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // For Triangle-to-Triangle distance, find minimum distance between all vertex pairs
        let vertices1 = [
            Point::new(self.first_coord().x(), self.first_coord().y()),
            Point::new(self.second_coord().x(), self.second_coord().y()),
            Point::new(self.third_coord().x(), self.third_coord().y()),
        ];

        let vertices2 = [
            Point::new(other.first_coord().x(), other.first_coord().y()),
            Point::new(other.second_coord().x(), other.second_coord().y()),
            Point::new(other.third_coord().x(), other.third_coord().y()),
        ];

        let mut min_distance = F::from(f64::INFINITY).unwrap();
        for &p1 in &vertices1 {
            for &p2 in &vertices2 {
                let dist = metric_space.distance(p1, p2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// Implementation for GeometryCollection with runtime type dispatch
impl<F, GC: GeometryCollectionTraitExt<T = F>> DistanceTrait<F, GeometryCollectionTag> for GC
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        // Find minimum distance between all geometry pairs from both collections
        let mut min_distance = F::from(f64::INFINITY).unwrap();

        for geom1 in self.geometries_ext() {
            for geom2 in other.geometries_ext() {
                let dist = geom1.distance_trait(metric_space, &geom2);
                if dist < min_distance {
                    min_distance = dist;
                }
            }
        }

        min_distance
    }
}

// Critical: GeometryTag implementation for WKB compatibility with runtime type dispatch
impl<F, G: GeometryTraitExt<T = F>> DistanceTrait<F, GeometryTag> for G
where
    F: CoordFloat,
{
    fn distance_trait(
        &self,
        metric_space: &impl Distance<F, Point<F>, Point<F>>,
        other: &Self,
    ) -> F {
        use geo_traits_ext::GeometryTypeExt;

        match (self.as_type_ext(), other.as_type_ext()) {
            (GeometryTypeExt::Point(p1), GeometryTypeExt::Point(p2)) => {
                p1.distance_trait(metric_space, p2)
            }
            (GeometryTypeExt::Line(l1), GeometryTypeExt::Line(l2)) => {
                l1.distance_trait(metric_space, l2)
            }
            (GeometryTypeExt::LineString(ls1), GeometryTypeExt::LineString(ls2)) => {
                ls1.distance_trait(metric_space, ls2)
            }
            (GeometryTypeExt::Polygon(poly1), GeometryTypeExt::Polygon(poly2)) => {
                poly1.distance_trait(metric_space, poly2)
            }
            (GeometryTypeExt::MultiPoint(mp1), GeometryTypeExt::MultiPoint(mp2)) => {
                mp1.distance_trait(metric_space, mp2)
            }
            (GeometryTypeExt::MultiLineString(mls1), GeometryTypeExt::MultiLineString(mls2)) => {
                mls1.distance_trait(metric_space, mls2)
            }
            (GeometryTypeExt::MultiPolygon(mpoly1), GeometryTypeExt::MultiPolygon(mpoly2)) => {
                mpoly1.distance_trait(metric_space, mpoly2)
            }
            (
                GeometryTypeExt::GeometryCollection(gc1),
                GeometryTypeExt::GeometryCollection(gc2),
            ) => gc1.distance_trait(metric_space, gc2),
            (GeometryTypeExt::Rect(r1), GeometryTypeExt::Rect(r2)) => {
                r1.distance_trait(metric_space, r2)
            }
            (GeometryTypeExt::Triangle(t1), GeometryTypeExt::Triangle(t2)) => {
                t1.distance_trait(metric_space, t2)
            }
            // For cross-type combinations, fall back to point-to-point distance
            _ => {
                // For simplicity, return zero for mixed types for now
                // TODO: Implement proper cross-type distance calculations
                F::zero()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coord, Euclidean, Line, LineString, Point, Polygon};

    #[test]
    fn point_distance_ext() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(3.0, 4.0);

        assert_eq!(p1.distance_ext(&Euclidean, &p2), 5.0);
    }

    #[test]
    fn linestring_distance_ext() {
        let ls1 = LineString::from(vec![coord!(x: 0.0, y: 0.0), coord!(x: 1.0, y: 0.0)]);
        let ls2 = LineString::from(vec![coord!(x: 2.0, y: 0.0), coord!(x: 3.0, y: 0.0)]);

        // Distance between closest points: (1.0, 0.0) to (2.0, 0.0) = 1.0
        assert_eq!(ls1.distance_ext(&Euclidean, &ls2), 1.0);
    }

    #[test]
    fn line_distance_ext() {
        let line1 = Line::new(coord!(x: 0.0, y: 0.0), coord!(x: 1.0, y: 0.0));
        let line2 = Line::new(coord!(x: 2.0, y: 0.0), coord!(x: 3.0, y: 0.0));

        // Distance between closest endpoints: (1.0, 0.0) to (2.0, 0.0) = 1.0
        assert_eq!(line1.distance_ext(&Euclidean, &line2), 1.0);
    }

    #[test]
    fn polygon_distance_ext() {
        let poly1 = Polygon::new(
            LineString::from(vec![
                coord!(x: 0.0, y: 0.0),
                coord!(x: 1.0, y: 0.0),
                coord!(x: 1.0, y: 1.0),
                coord!(x: 0.0, y: 1.0),
                coord!(x: 0.0, y: 0.0),
            ]),
            vec![],
        );
        let poly2 = Polygon::new(
            LineString::from(vec![
                coord!(x: 2.0, y: 0.0),
                coord!(x: 3.0, y: 0.0),
                coord!(x: 3.0, y: 1.0),
                coord!(x: 2.0, y: 1.0),
                coord!(x: 2.0, y: 0.0),
            ]),
            vec![],
        );

        // Distance between closest vertices: (1.0, 0.0) to (2.0, 0.0) = 1.0
        assert_eq!(poly1.distance_ext(&Euclidean, &poly2), 1.0);
    }
}
