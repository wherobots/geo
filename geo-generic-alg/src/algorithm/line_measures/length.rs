use super::Distance;
use crate::{CoordFloat, Line, LineString, MultiLineString, Point};
use geo_traits::{CoordTrait, PolygonTrait};
use geo_traits_ext::*;

/// Calculate the length of a geometry using a given [metric space](crate::algorithm::line_measures::metric_spaces).
///
/// For 1D geometries (Line, LineString, MultiLineString), this returns the actual length.
/// For 2D geometries (Polygon, MultiPolygon, Rect, Triangle), this returns the perimeter.
/// For 0D geometries (Point, MultiPoint), this returns zero.
///
/// # Examples
/// ```
/// use geo::algorithm::line_measures::{Length, Euclidean, Haversine};
///
/// let line_string = geo::wkt!(LINESTRING(
///     0.0 0.0,
///     3.0 4.0,
///     3.0 5.0
/// ));
/// assert_eq!(Euclidean.length(&line_string), 6.);
///
/// let line_string_lon_lat = geo::wkt!(LINESTRING (
///     -47.9292 -15.7801f64,
///     -58.4173 -34.6118,
///     -70.6483 -33.4489
/// ));
/// assert_eq!(Haversine.length(&line_string_lon_lat).round(), 3_474_956.0);
/// ```
pub trait Length<F: CoordFloat> {
    fn length(&self, geometry: &impl LengthMeasurable<F>) -> F;
}

/// Something which can be measured by a [metric space](crate::algorithm::line_measures::metric_spaces).
///
/// For 1D geometries, returns the length. For 2D geometries, returns the perimeter.
///
/// It's typically more convenient to use the [`Length`] trait instead of this trait directly.
///
/// # Examples
/// ```
/// use geo::algorithm::line_measures::{LengthMeasurable, Euclidean, Haversine};
///
/// let line_string = geo::wkt!(LINESTRING(
///     0.0 0.0,
///     3.0 4.0,
///     3.0 5.0
/// ));
/// assert_eq!(line_string.length(&Euclidean), 6.);
///
/// let line_string_lon_lat = geo::wkt!(LINESTRING (
///     -47.9292 -15.7801f64,
///     -58.4173 -34.6118,
///     -70.6483 -33.4489
/// ));
/// assert_eq!(line_string_lon_lat.length(&Haversine).round(), 3_474_956.0);
/// ```
pub trait LengthMeasurable<F: CoordFloat> {
    fn length(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F;
}

impl<F: CoordFloat, PointDistance: Distance<F, Point<F>, Point<F>>> Length<F> for PointDistance {
    fn length(&self, geometry: &impl LengthMeasurable<F>) -> F {
        geometry.length(self)
    }
}

impl<F: CoordFloat> LengthMeasurable<F> for Line<F> {
    fn length(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        metric_space.distance(self.start_point(), self.end_point())
    }
}

impl<F: CoordFloat> LengthMeasurable<F> for LineString<F> {
    fn length(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        let mut length = F::zero();
        for line in self.lines() {
            length = length + line.length(metric_space);
        }
        length
    }
}

impl<F: CoordFloat> LengthMeasurable<F> for MultiLineString<F> {
    fn length(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        let mut length = F::zero();
        for line in self {
            length = length + line.length(metric_space);
        }
        length
    }
}

/// Extension trait that enables the modern Length API for WKB and other generic geometry types.
///
/// This provides the same API as the concrete `LengthMeasurable` implementations but works with
/// any geometry type that implements the geo-traits-ext pattern.
///
/// For 1D geometries (Line, LineString, MultiLineString), returns the length.
/// For 2D geometries (Polygon, MultiPolygon, Rect, Triangle), returns the perimeter.
/// For 0D geometries (Point, MultiPoint), returns zero.
///
/// # Examples
/// ```
/// use geo_generic_alg::algorithm::line_measures::{LengthMeasurableExt, Euclidean};
/// // Example WKB bytes for LINESTRING(0 0, 3 4, 3 5) in little-endian
/// let wkb_bytes: &[u8] = &[
///     1, // little endian
///     2, 0, 0, 0, // geometry type: LineString (2)
///     3, 0, 0, 0, // num points: 3
///     0, 0, 0, 0, 0, 0, 0, 0, // x0 = 0.0
///     0, 0, 0, 0, 0, 0, 0, 0, // y0 = 0.0
///     0, 0, 0, 0, 0, 0, 8, 64, // x1 = 3.0
///     0, 0, 0, 0, 0, 0, 16, 64, // y1 = 4.0
///     0, 0, 0, 0, 0, 0, 8, 64, // x2 = 3.0
///     0, 0, 0, 0, 0, 0, 20, 64, // y2 = 5.0
/// ];
/// let wkb_geom = geo_generic_tests::wkb::reader::read_wkb(wkb_bytes).unwrap();
/// let length = wkb_geom.length_ext(&Euclidean);
/// assert_eq!(length, 6.0);
/// ```
pub trait LengthMeasurableExt<F: CoordFloat> {
    /// Calculate the length using the given metric space.
    fn length_ext(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F;
}

// Implementation for WKB and other generic geometries using the type-tag pattern
impl<F, G> LengthMeasurableExt<F> for G
where
    F: CoordFloat,
    G: GeoTraitExtWithTypeTag + LengthMeasurableTrait<F, G::Tag>,
{
    fn length_ext(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        self.length_trait(metric_space)
    }
}

// Internal trait that handles the actual length computation for different geometry types
trait LengthMeasurableTrait<F, GT: GeoTypeTag>
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F;
}

// Implementation for Line geometries
impl<F, L: LineTraitExt<T = F>> LengthMeasurableTrait<F, LineTag> for L
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        let start = Point::new(self.start_coord().x, self.start_coord().y);
        let end = Point::new(self.end_coord().x, self.end_coord().y);
        metric_space.distance(start, end)
    }
}

// Implementation for LineString geometries
impl<F, LS: LineStringTraitExt<T = F>> LengthMeasurableTrait<F, LineStringTag> for LS
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        let mut length = F::zero();
        for line in self.lines() {
            let start = Point::new(line.start_coord().x, line.start_coord().y);
            let end = Point::new(line.end_coord().x, line.end_coord().y);
            length = length + metric_space.distance(start, end);
        }
        length
    }
}

// Implementation for MultiLineString geometries
impl<F, MLS: MultiLineStringTraitExt<T = F>> LengthMeasurableTrait<F, MultiLineStringTag> for MLS
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        let mut length = F::zero();
        for line_string in self.line_strings_ext() {
            length = length + line_string.length_trait(metric_space);
        }
        length
    }
}

// For geometry types that don't have a meaningful length (return zero)
impl<F, P: PointTraitExt<T = F>> LengthMeasurableTrait<F, PointTag> for P
where
    F: CoordFloat,
{
    fn length_trait(&self, _metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        F::zero()
    }
}

impl<F, MP: MultiPointTraitExt<T = F>> LengthMeasurableTrait<F, MultiPointTag> for MP
where
    F: CoordFloat,
{
    fn length_trait(&self, _metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        F::zero()
    }
}

// Helper function to calculate the perimeter of a linestring using a metric space
fn linestring_perimeter_with_metric<F, LS: LineStringTraitExt<T = F>>(
    linestring: &LS,
    metric_space: &impl Distance<F, Point<F>, Point<F>>,
) -> F
where
    F: CoordFloat,
{
    let mut perimeter = F::zero();
    for line in linestring.lines() {
        let start_coord = line.start_coord();
        let end_coord = line.end_coord();
        let start_point = Point::new(start_coord.x(), start_coord.y());
        let end_point = Point::new(end_coord.x(), end_coord.y());
        perimeter = perimeter + metric_space.distance(start_point, end_point);
    }
    perimeter
}

// Helper function to calculate the perimeter of a ring using the basic LineStringTrait
fn ring_perimeter_with_metric<F, LS>(
    ring: &LS,
    metric_space: &impl Distance<F, Point<F>, Point<F>>,
) -> F
where
    F: CoordFloat,
    LS: geo_traits::LineStringTrait<T = F>,
{
    let mut perimeter = F::zero();
    let num_coords = ring.num_coords();
    if num_coords > 1 {
        for i in 0..(num_coords - 1) {
            let start_coord = ring.coord(i).unwrap();
            let end_coord = ring.coord(i + 1).unwrap();
            let start_point = Point::new(start_coord.x(), start_coord.y());
            let end_point = Point::new(end_coord.x(), end_coord.y());
            perimeter = perimeter + metric_space.distance(start_point, end_point);
        }
    }
    perimeter
}

impl<F, P: PolygonTraitExt<T = F>> LengthMeasurableTrait<F, PolygonTag> for P
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        // For polygons, return the perimeter (length of the boundary)
        let mut total_perimeter = match self.exterior_ext() {
            Some(exterior) => linestring_perimeter_with_metric(&exterior, metric_space),
            None => F::zero(),
        };

        // Add interior rings perimeter
        for interior in self.interiors_ext() {
            total_perimeter =
                total_perimeter + linestring_perimeter_with_metric(&interior, metric_space);
        }

        total_perimeter
    }
}

impl<F, MP: MultiPolygonTraitExt<T = F>> LengthMeasurableTrait<F, MultiPolygonTag> for MP
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        // For multipolygons, return the sum of all polygon perimeters
        let mut total_perimeter = F::zero();
        for polygon in self.polygons() {
            // Calculate perimeter for each polygon
            let mut polygon_perimeter = match polygon.exterior() {
                Some(exterior) => ring_perimeter_with_metric(&exterior, metric_space),
                None => F::zero(),
            };

            // Add interior rings perimeter
            for interior in polygon.interiors() {
                polygon_perimeter =
                    polygon_perimeter + ring_perimeter_with_metric(&interior, metric_space);
            }

            total_perimeter = total_perimeter + polygon_perimeter;
        }
        total_perimeter
    }
}

impl<F, R: RectTraitExt<T = F>> LengthMeasurableTrait<F, RectTag> for R
where
    F: CoordFloat,
{
    fn length_trait(&self, _metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        // For rectangles, return the perimeter (use width/height methods if available)
        let width = self.width();
        let height = self.height();
        let two = F::one() + F::one();
        two * (width + height)
    }
}

impl<F, T: TriangleTraitExt<T = F>> LengthMeasurableTrait<F, TriangleTag> for T
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        // For triangles, return the perimeter (sum of all three sides)
        let coord0 = self.first_coord();
        let coord1 = self.second_coord();
        let coord2 = self.third_coord();

        let p0 = Point::new(coord0.x(), coord0.y());
        let p1 = Point::new(coord1.x(), coord1.y());
        let p2 = Point::new(coord2.x(), coord2.y());

        let side1 = metric_space.distance(p0, p1);
        let side2 = metric_space.distance(p1, p2);
        let side3 = metric_space.distance(p2, p0);

        side1 + side2 + side3
    }
}

// Implementation for GeometryCollection with runtime type dispatch
impl<F, GC: GeometryCollectionTraitExt<T = F>> LengthMeasurableTrait<F, GeometryCollectionTag>
    for GC
where
    F: CoordFloat,
{
    fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F {
        self.geometries_ext()
            .map(|g| match g.as_type_ext() {
                GeometryTypeExt::Point(_) => F::zero(),
                GeometryTypeExt::Line(line) => line.length_trait(metric_space),
                GeometryTypeExt::LineString(ls) => ls.length_trait(metric_space),
                GeometryTypeExt::Polygon(polygon) => polygon.length_trait(metric_space),
                GeometryTypeExt::MultiPoint(_) => F::zero(),
                GeometryTypeExt::MultiLineString(mls) => mls.length_trait(metric_space),
                GeometryTypeExt::MultiPolygon(mp) => mp.length_trait(metric_space),
                GeometryTypeExt::GeometryCollection(gc) => gc.length_trait(metric_space),
                GeometryTypeExt::Rect(rect) => rect.length_trait(metric_space),
                GeometryTypeExt::Triangle(triangle) => triangle.length_trait(metric_space),
            })
            .fold(F::zero(), |acc, next| acc + next)
    }
}

// Critical: GeometryTag implementation for WKB compatibility
impl<F, G: GeometryTraitExt<T = F>> LengthMeasurableTrait<F, GeometryTag> for G
where
    F: CoordFloat,
{
    // This macro delegates the `length_trait` method to the appropriate geometry variant.
    // It is critical for WKB (Well-Known Binary) compatibility, ensuring that trait methods
    // are correctly dispatched for all geometry types when deserializing from WKB.
    crate::geometry_trait_ext_delegate_impl! {
        fn length_trait(&self, metric_space: &impl Distance<F, Point<F>, Point<F>>) -> F;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{coord, Euclidean, Geodesic, Haversine, Rhumb};

    #[test]
    fn lines() {
        // london to paris
        let line = Line::new(
            coord!(x: -0.1278f64, y: 51.5074),
            coord!(x: 2.3522, y: 48.8566),
        );

        assert_eq!(
            343_923., // meters
            Geodesic.length(&line).round()
        );
        assert_eq!(
            343_572., // meters
            Rhumb.length(&line).round()
        );
        assert_eq!(
            343_557., // meters
            Haversine.length(&line).round()
        );

        // computing Euclidean length of an unprojected (lng/lat) line gives a nonsense answer
        assert_eq!(
            4., // nonsense!
            Euclidean.length(&line).round()
        );
        // london to paris in EPSG:3035
        let projected_line = Line::new(
            coord!(x: 3620451.74f64, y: 3203901.44),
            coord!(x: 3760771.86, y: 2889484.80),
        );
        assert_eq!(344_307., Euclidean.length(&projected_line).round());
    }

    #[test]
    fn line_strings() {
        let line_string = LineString::new(vec![
            coord!(x: -58.3816f64, y: -34.6037), // Buenos Aires, Argentina
            coord!(x: -77.0428, y: -12.0464),    // Lima, Peru
            coord!(x: -47.9292, y: -15.7801),    // Brasília, Brazil
        ]);

        assert_eq!(
            6_302_220., // meters
            Geodesic.length(&line_string).round()
        );
        assert_eq!(
            6_308_683., // meters
            Rhumb.length(&line_string).round()
        );
        assert_eq!(
            6_304_387., // meters
            Haversine.length(&line_string).round()
        );

        // computing Euclidean length of an unprojected (lng/lat) gives a nonsense answer
        assert_eq!(
            59., // nonsense!
            Euclidean.length(&line_string).round()
        );
        // EPSG:102033
        let projected_line_string = LineString::from(vec![
            coord!(x: 143042.46f64, y: -1932485.45), // Buenos Aires, Argentina
            coord!(x: -1797084.08, y: 583528.84),    // Lima, Peru
            coord!(x: 1240052.27, y: 207169.12),     // Brasília, Brazil
        ]);
        assert_eq!(6_237_538., Euclidean.length(&projected_line_string).round());
    }

    // Tests for LengthMeasurableExt - adapted from euclidean_length.rs
    mod length_measurable_ext_tests {
        use super::*;
        use crate::{
            coord, line_string, polygon, Geometry, GeometryCollection, Line, MultiLineString,
            MultiPoint, MultiPolygon, Point, Polygon,
        };

        #[test]
        fn empty_linestring_test() {
            let linestring = line_string![];
            assert_relative_eq!(0.0_f64, linestring.length_ext(&Euclidean));
        }

        #[test]
        fn linestring_one_point_test() {
            let linestring = line_string![(x: 0., y: 0.)];
            assert_relative_eq!(0.0_f64, linestring.length_ext(&Euclidean));
        }

        #[test]
        fn linestring_test() {
            let linestring = line_string![
                (x: 1., y: 1.),
                (x: 7., y: 1.),
                (x: 8., y: 1.),
                (x: 9., y: 1.),
                (x: 10., y: 1.),
                (x: 11., y: 1.)
            ];
            assert_relative_eq!(10.0_f64, linestring.length_ext(&Euclidean));
        }

        #[test]
        fn multilinestring_test() {
            let mline = MultiLineString::new(vec![
                line_string![
                    (x: 1., y: 0.),
                    (x: 7., y: 0.),
                    (x: 8., y: 0.),
                    (x: 9., y: 0.),
                    (x: 10., y: 0.),
                    (x: 11., y: 0.)
                ],
                line_string![
                    (x: 0., y: 0.),
                    (x: 0., y: 5.)
                ],
            ]);
            assert_relative_eq!(15.0_f64, mline.length_ext(&Euclidean));
        }

        #[test]
        fn line_test() {
            let line0 = Line::new(coord! { x: 0., y: 0. }, coord! { x: 0., y: 1. });
            let line1 = Line::new(coord! { x: 0., y: 0. }, coord! { x: 3., y: 4. });
            assert_relative_eq!(line0.length_ext(&Euclidean), 1.);
            assert_relative_eq!(line1.length_ext(&Euclidean), 5.);
        }

        #[test]
        fn polygon_perimeter_test() {
            let polygon: Polygon<f64> = polygon![
                (x: 0., y: 0.),
                (x: 4., y: 0.),
                (x: 4., y: 4.),
                (x: 0., y: 4.),
                (x: 0., y: 0.),
            ];
            // For polygons, length_ext returns the perimeter: 4 + 4 + 4 + 4 = 16
            assert_relative_eq!(polygon.length_ext(&Euclidean), 16.0);
        }

        #[test]
        fn point_returns_zero_test() {
            let point = Point::new(3.0, 4.0);
            // Points have no length dimension
            assert_relative_eq!(point.length_ext(&Euclidean), 0.0);
        }

        #[test]
        fn comprehensive_test_scenarios() {
            // Test cases matching the Python pytest scenarios

            // LINESTRING EMPTY
            let empty_linestring: crate::LineString<f64> = line_string![];
            assert_relative_eq!(empty_linestring.length_ext(&Euclidean), 0.0);

            // POINT (0 0)
            let point = Point::new(0.0, 0.0);
            assert_relative_eq!(point.length_ext(&Euclidean), 0.0);

            // LINESTRING (0 0, 0 1) - length should be 1
            let linestring = line_string![(x: 0., y: 0.), (x: 0., y: 1.)];
            assert_relative_eq!(linestring.length_ext(&Euclidean), 1.0);

            // MULTIPOINT ((0 0), (1 1)) - should be 0
            let multipoint = MultiPoint::new(vec![Point::new(0.0, 0.0), Point::new(1.0, 1.0)]);
            assert_relative_eq!(multipoint.length_ext(&Euclidean), 0.0);

            // MULTILINESTRING ((0 0, 1 1), (1 1, 2 2)) - should be ~2.828427
            let multilinestring = MultiLineString::new(vec![
                line_string![(x: 0., y: 0.), (x: 1., y: 1.)],
                line_string![(x: 1., y: 1.), (x: 2., y: 2.)],
            ]);
            assert_relative_eq!(
                multilinestring.length_ext(&Euclidean),
                2.8284271247461903,
                epsilon = 1e-10
            );

            // POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0)) - should be 4.0 (perimeter: 1+1+1+1)
            let polygon = polygon![
                (x: 0., y: 0.),
                (x: 1., y: 0.),
                (x: 1., y: 1.),
                (x: 0., y: 1.),
                (x: 0., y: 0.),
            ];
            assert_relative_eq!(polygon.length_ext(&Euclidean), 4.0);

            // MULTIPOLYGON - should be 8.0 (two polygons with perimeter 4.0 each)
            let multipolygon = MultiPolygon::new(vec![
                polygon![
                    (x: 0., y: 0.),
                    (x: 1., y: 0.),
                    (x: 1., y: 1.),
                    (x: 0., y: 1.),
                    (x: 0., y: 0.),
                ],
                polygon![
                    (x: 2., y: 2.),
                    (x: 3., y: 2.),
                    (x: 3., y: 3.),
                    (x: 2., y: 3.),
                    (x: 2., y: 2.),
                ],
            ]);
            assert_relative_eq!(multipolygon.length_ext(&Euclidean), 8.0);

            // GEOMETRYCOLLECTION (LINESTRING (0 0, 1 1), POLYGON (...), LINESTRING (0 0, 1 1))
            // Should sum linestrings + polygon perimeter: 2 * sqrt(2) + 4 ≈ 6.8284271247461903
            let collection = GeometryCollection::new_from(vec![
                Geometry::LineString(line_string![(x: 0., y: 0.), (x: 1., y: 1.)]), // sqrt(2) ≈ 1.414
                Geometry::Polygon(polygon![
                    (x: 0., y: 0.),
                    (x: 1., y: 0.),
                    (x: 1., y: 1.),
                    (x: 0., y: 1.),
                    (x: 0., y: 0.),
                ]), // perimeter = 4.0
                Geometry::LineString(line_string![(x: 0., y: 0.), (x: 1., y: 1.)]), // sqrt(2) ≈ 1.414
            ]);
            assert_relative_eq!(
                collection.length_ext(&Euclidean),
                2.8284271247461903 + 4.0, // 2*sqrt(2) + 4
                epsilon = 1e-10
            );
        }

        #[test]
        fn test_different_metric_spaces() {
            // Test that different metric spaces work with LengthMeasurableExt
            let linestring = line_string![(x: 0., y: 0.), (x: 3., y: 4.)];

            // Euclidean should give us 5.0
            assert_relative_eq!(linestring.length_ext(&Euclidean), 5.0);

            // Test with geographic coordinates (lon/lat) - these will give nonsense results
            // but should demonstrate the API works with different metric spaces
            let lon_lat_line = line_string![
                (x: -0.1278f64, y: 51.5074), // London
                (x: 2.3522, y: 48.8566)      // Paris
            ];

            // These should all work without errors (values are from the existing tests)
            assert_eq!(Haversine.length(&lon_lat_line).round(), 343_557.);
            assert_eq!(lon_lat_line.length_ext(&Haversine).round(), 343_557.);

            // Verify both APIs give the same result
            assert_relative_eq!(
                Haversine.length(&lon_lat_line),
                lon_lat_line.length_ext(&Haversine),
                epsilon = 1e-6
            );
        }

        #[test]
        fn test_polygon_with_holes() {
            // Test polygon with interior rings (holes)
            let polygon = Polygon::new(
                LineString::new(vec![
                    coord! { x: 0., y: 0. },
                    coord! { x: 10., y: 0. },
                    coord! { x: 10., y: 10. },
                    coord! { x: 0., y: 10. },
                    coord! { x: 0., y: 0. },
                ]),
                vec![LineString::new(vec![
                    coord! { x: 2., y: 2. },
                    coord! { x: 8., y: 2. },
                    coord! { x: 8., y: 8. },
                    coord! { x: 2., y: 8. },
                    coord! { x: 2., y: 2. },
                ])],
            );
            // Exterior perimeter: 40 (10+10+10+10), Interior perimeter: 24 (6+6+6+6)
            assert_relative_eq!(polygon.length_ext(&Euclidean), 64.0);
        }

        #[test]
        fn test_triangle_perimeter() {
            use crate::Triangle;
            // Right triangle with sides 3, 4, 5
            let triangle = Triangle::new(
                coord! { x: 0., y: 0. },
                coord! { x: 3., y: 0. },
                coord! { x: 0., y: 4. },
            );
            // Perimeter should be 3 + 4 + 5 = 12
            assert_relative_eq!(triangle.length_ext(&Euclidean), 12.0);
        }

        #[test]
        fn test_rect_perimeter() {
            use crate::Rect;
            // Rectangle 3x4
            let rect = Rect::new(coord! { x: 0., y: 0. }, coord! { x: 3., y: 4. });
            // Perimeter should be 2*(3+4) = 14
            assert_relative_eq!(rect.length_ext(&Euclidean), 14.0);
        }
    }
}
