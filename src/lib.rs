#[macro_use]
extern crate nom;
extern crate geo;
extern crate dbf;

use geo::*;

use nom::*;
use dbf::DbfFile;
use dbf::Record as DbfRecord;

use std::path::Path;
use std::io::prelude::*;
use std::io::SeekFrom;
use std::fs::File;

#[derive(Debug)]
struct Header {
    version: i32,
    shape_type: i32,
    mbr: [f64; 4],
    zrange: [f64; 2],
    mrange: [f64; 2],
}

fn read_bytes<R: Read+Seek>(input: &mut R, start: u64, length: usize) -> Result<Vec<u8>, String> {
    //println!("Want to read from start={} for length={} bytes to end={}", start, length, start+length as u64);
    let mut res_vec = vec![0; length];

    try!(input.seek(SeekFrom::Start(start)).map_err(|_| "couldn't seek".to_string()));

    try!(input.read_exact(&mut res_vec).map_err(|_| "Couldn't read bytes".to_string()));
    Ok(res_vec)
}

named!(parse_header<Header>,
   do_parse!(
       // TODO replace this with the tag 9994 (\x00\x00\x27\x0a)?
       take!(4) >>
       take!(5*4)       >>
       file_length: take!(4)       >>
       version: le_i32 >>
       shape_type: le_i32 >>
       mbr: count_fixed!( f64, le_f64, 4) >>
       zrange: count_fixed!( f64, le_f64, 2) >>
       mrange: count_fixed!( f64, le_f64, 2) >>

       ( Header{ version: version, shape_type: shape_type, mbr: mbr, zrange: zrange, mrange: mrange } )

   )
);

named!(parse_coordinate<Coordinate<f64> >,
       do_parse!(
           x: le_f64 >>
           y: le_f64 >>
           ( Coordinate{ x: x, y: y} )
       )
);

named!(parse_point<Point<f64> >,
       do_parse!(
           x: le_f64 >> y: le_f64 >>
           ( Point::new(x, y) )
       )
);

named!(parse_bbox<Bbox<f64> >,
       do_parse!(
           xmin: le_f64 >> ymin: le_f64 >> xmax: le_f64 >> ymax: le_f64 >>
           ( Bbox{ xmin: xmin, ymin: ymin, xmax: xmax, ymax: ymax } )
       )
      );

named!(parse_shx_offset<(i32, i32)>, tuple!( be_i32, be_i32 ) );

/// Returns true iff the ring is a clockwise ring. false if anti-clockwise. presumes a valid ring
fn ring_is_clockwise(ring: &LineString<f64>) -> bool {
    // FIXME can we use .windows here instead?
    let res: f64 = (0..ring.0.len()-1).map(|i| ( ring.0[i].x()*ring.0[i+1].y()  - ring.0[i+1].x()*ring.0[i].y() ) ).sum();
    res < 0.
}

fn parse_polygon(i: &[u8]) -> IResult<&[u8], Geometry<f64>> {
    let (i, _bbox) = try_parse!(i, parse_bbox);
    let (i, num_parts) = try_parse!(i, le_i32);
    let (i, num_points) = try_parse!(i, le_i32);
    let (i, mut parts) = try_parse!(i, count!( le_i32, num_parts as usize ) );
    parts.push(num_points);
    let parts: Vec<usize> = parts.into_iter().map(|x| x as usize).collect();
    //println!("Parts: {:?}", parts);

    let (i, points) = try_parse!(i, count!( parse_point, num_points as usize ) );
    let points = points.as_slice();
    //println!("Bbox {:?} num_parts {} num_points {} parts {:?} points {:?}", _bbox, num_parts, num_points, parts, points);

    let mut polygons: Vec<Polygon<f64>> = Vec::new();
    for ring_id in 0..num_parts {
        let ring_id = ring_id as usize;
        let mut ring = Vec::new();   // TODO with_capacity
        ring.extend_from_slice(&points[ parts[ring_id] .. parts[ring_id+1]-1 ]);
        let linestring = LineString(ring);
        if ring_is_clockwise(&linestring) {
            // new polygon
            polygons.push(Polygon::new(linestring, vec![]));
        } else {
            // this is an interior
            if polygons.len() > 0 {
                let last_idx = polygons.len()-1;
                polygons[last_idx].interiors.push(linestring);
            } else {
                // FIXME
                // TODO Should we just presume that this is an exterious ring, even if it's
                // counter-clockwise?
                //println!("Have 'interior' polygon with no previous exterior ring??");
                //println!("Bbox {:?} num_parts {} num_points {} parts {:?} points {:?}", _bbox, num_parts, num_points, parts, points);
                //panic!();
            }
        }
    }

    let result: Geometry<f64> = if polygons.len() == 1 {
        // FIXME Is swap_remove the best here? I want to remove object at 0. I don't need polygons
        // afterwards
        Geometry::Polygon(polygons.swap_remove(0))
    } else {
        Geometry::MultiPolygon(MultiPolygon(polygons))
    };

    IResult::Done(i, result)
}


fn parse_record(i: &[u8]) -> IResult<&[u8], (i32, Option<Geometry<f64>>)> {
    // header
    let (i, rec_num) = try_parse!(i, be_i32);
    let (i, _rec_len) = try_parse!(i, be_i32);
    let (i, shape_type) = try_parse!(i, le_i32);

    // FIXME there is a bug here, shape_bytes isn't big enough
    //let (i, shape_bytes) = try_parse!(i, take!(2*(rec_len-4)));
 
    let (i, shape) = match shape_type {
        0 => { (i, None) },
        1 => { parse_point(&i).map(|p| Some(Geometry::Point(p))).unwrap() },
        //3 => { parse_polyline(&shape_bytes).unwrap() },
        5 => { parse_polygon(&i).map(|p| Some(p)).unwrap() },
        //8 => { parse_multipoint(&shape_bytes).unwrap() },
        _ => { (i, None) },
    };

    IResult::Done(&i, (rec_num, shape))
}

named!(parse_records<Vec<(i32, Option<Geometry<f64> >) > >, many0!( parse_record ) );


pub struct Shapefile {
    _shp_file_handle: File,
    _shx_file_handle: File,

    _dbf_file: DbfFile<File>,

    _header: Header,
}

#[derive(Debug)]
pub struct Record {
    id: i32,
    geometry: Option<Geometry<f64>>,
    attributes: DbfRecord,
}

impl Shapefile {
    pub fn open(filename: &Path) -> Self {
        let mut shp_file = File::open(filename).unwrap();
        let shx_filename = filename.with_extension("shx");

        let shx_file = File::open(shx_filename).unwrap();

        let header_bytes = read_bytes(&mut shp_file, 0, 100).unwrap();

        let hdr: Header = parse_header(&header_bytes).to_result().unwrap();

        let dbf_file = DbfFile::open_file(&filename.with_extension("dbf")).unwrap();

        Shapefile{ _header: hdr, _shp_file_handle: shp_file, _shx_file_handle: shx_file, _dbf_file: dbf_file }
    }

    pub fn read_all(&mut self) -> Vec<Record> {
        self._shp_file_handle.seek(SeekFrom::Start(0)).unwrap();
        let mut bytes = Vec::new();
        self._shp_file_handle.read_to_end(&mut bytes).unwrap();

        let (bytes, hdr) = parse_header(&bytes).unwrap();
        println!("Header {:?}", hdr);

        let (bytes, recs) = parse_records(&bytes).unwrap();
        let mut objects = Vec::with_capacity(recs.len());
        for (idx, rec) in recs.into_iter().enumerate() {
            let (id, geometry) = rec;
            //println!("idx {} id {}", id
            let attributes = self._dbf_file.record(idx as u32).unwrap();
            objects.push(Record{ id: id, geometry: geometry, attributes: attributes });
        }

        objects
    }

    pub fn record(&mut self, rec_id: u32) -> Option<Record> {
        let offset = (100 + 8 * rec_id) as u64;
        let shx_bytes = read_bytes(&mut self._shx_file_handle, offset, 8).unwrap();
        let shp_start_offset = parse_shx_offset(&shx_bytes).to_result().unwrap();
        let shp_bytes = read_bytes(&mut self._shp_file_handle, (shp_start_offset.0) as u64, (shp_start_offset.1) as usize).unwrap();

        let (id, geometry): (i32, Option<Geometry<f64>>) = parse_record(&shp_bytes).to_result().unwrap();
        let attributes = self._dbf_file.record(rec_id).unwrap();
        let shp_record = Record{ id: id, geometry: geometry, attributes: attributes };

        Some(shp_record)

    }

    pub fn records(self) -> ShapefileIterator {
        ShapefileIterator{ _shapefile: self, _next_rec: 0, _done: false }
    }
}

pub struct ShapefileIterator {
    _shapefile: Shapefile,
    _next_rec: u32,
    _done: bool,
}

impl ShapefileIterator {
    pub fn into_inner(self) -> Shapefile {
        self._shapefile
    }
}

impl Iterator for ShapefileIterator {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        if self._done {
            None
        } else {
            let rec = self._shapefile.record(self._next_rec);
            self._next_rec = self._next_rec + 1;
            rec
        }
    }
}
