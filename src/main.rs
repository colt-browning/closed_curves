/*
Ported to Rust and extended with new functionality in 2025 by Andrei Zabolotskii.
Original notice:
   Program by F.S.Duzhin, fduzhin@mccme.ru, January 1998.
   Ported to GNU Pascal by S.V.Duzhin.

   Enumeration of closed curves in the plane in the sense of 
   V.I.Arnold "Topological Invariants of Plane Curves and Caustics", 
   University Lecture Series, v.5, AMS 1994.

   Algorithm of counting long curves by S.M.Gusein-Zade (1994)
   Algorithm of counting closed curves by F.S.Duzhin (1998)

   Published in:
   S.M.Gusein-Zade, F.S.Duzhin. On the number of topological types 
   of plane curves. (Russian) Uspekhi Mat. Nauk 53 (1998),
   no. 3(321), 197--198. 
*/

const MAX_X: usize = 12; // can be set to 15

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Arc {
	w: usize,
	up: bool,
}

#[derive(Debug, Clone, Default)]
struct Curve {
	nx: usize,
	arcs: [Arc; 2*MAX_X + 1],
	dir_r: bool,
}

impl PartialEq for Curve {
	fn eq(&self, other: &Curve) -> bool {
		self.nx == other.nx && self.arcs[..=2*self.nx] == other.arcs[..=2*other.nx]
	}
}

impl Eq for Curve {}

use std::fmt;

impl fmt::Display for Curve {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "{}[ ", if self.dir_r {"→"} else {"←"})?;
		for i in 1..=(2*self.nx) {
			write!(f, "{}{} ", self.arcs[i].w, if self.arcs[i].up {"↑"} else {"↓"})?;
		}
		write!(f, "]")
	}
}

type Numeration = [usize; 2*MAX_X + 2];
type Free = [bool; 2*MAX_X + 2];

impl Curve {
	fn _check(&self) {
		for i in 1..=(2*self.nx) {
			let a1 = self.arcs[i];
			let a2 = self.arcs[a1.w];
			assert_eq!(a2.w, i, "{self}");
			assert_eq!(a1.up, a2.up, "{self}");
		}
	}
	
	fn _check_parity(&self) -> bool {
		let g = self._gauss_paired();
		for i in 1..=(2*self.nx) {
			let j = g[i];
			if (j as isize - i as isize) % 2 == 0 {
				return false
			}
		}
		true
	}
	
	fn put_indices(&self, n: &mut Numeration) {
		let mut i = 2*self.nx;
		let mut k = 1;
		n[i+1] = k;
		while i > 0 {
			k += 1;
			i = self.arcs[i].w + 1;
			n[i] = k;
			k += 1;
			i = self.arcs[i].w;
			n[i] = k;
			i -= 1;
		}
	}
	
	fn look_for_free_segments(&self, connect: &mut Free, side: bool) {
		let (mut l, mut r) = (1, 1);
		connect[1] = true;
		connect[2*self.nx+1] = true;
		while r <= 2*self.nx {
			if self.arcs[r].up == side && self.arcs[r].w > r {
				r = self.arcs[r].w;
				while l <= r {
					l += 1;
					connect[l] = false;
				}
				l = r;
			} else {
				r += 1;
				l += 1;
				connect[l] = true;
			}
		}
	}
	
	fn make_arc(&mut self, h: &Curve, segment: usize, side: bool) {
		self.nx = h.nx + 1;
		for i in ((segment+1)..(2*self.nx)).rev() {
			self.arcs[i].w = h.arcs[i-1].w + if h.arcs[i-1].w >= segment { 1 } else { 0 };
			self.arcs[i].up = h.arcs[i-1].up;
		}
		for i in 1..segment {
			self.arcs[i].w = h.arcs[i].w + if h.arcs[i].w >= segment { 1 } else { 0 };
			self.arcs[i].up = h.arcs[i].up;
		}
		self.arcs[segment].w = 2*self.nx;
		self.arcs[segment].up = side;
		self.arcs[2*self.nx].w = segment;
		self.arcs[2*self.nx].up = side;
	}
	
	fn remake_curve(&mut self, ind: &Numeration) {
		let ha = self.arcs.clone();
		let cut = if self.arcs[2*self.nx].up {
			2 * self.nx - 1
		} else {
			self.arcs[2*self.nx].w - 1
		};
		if ind[cut+1]%2 == 0 {
			self.dir_r = !self.dir_r;
		}
		for i in 1..=(self.nx*2) {
			let t = &mut self.arcs[if i > cut { i-cut } else { i+2*self.nx-cut }];
			if ha[i].w <= cut {
				t.w = ha[i].w+2*self.nx-cut;
			} else {
				t.w = ha[i].w-cut;
			}
			t.up = ha[i].up;
		}
	}
	
	fn kill_arc(&mut self, ind: &Numeration) -> (usize, bool) {
		let h = self.clone();
		if ind[self.arcs[2*self.nx].w]%2 == 1 {
			self.nx -= 1;
			for i in 1..h.arcs[2*h.nx].w {
				self.arcs[i] = h.arcs[i];
				if h.arcs[i].w >= h.arcs[2*h.nx].w { self.arcs[i].w -= 1 };
			}
			for i in (h.arcs[2*h.nx].w+1)..(2*h.nx) {
				self.arcs[i-1] = h.arcs[i];
				if h.arcs[i].w >= h.arcs[2*h.nx].w { self.arcs[i-1].w -= 1 };
			}
			(ind[2*h.nx] - 1, h.arcs[2*h.nx].up)
		} else {
			let border = h.arcs[2*h.nx].w;
			let side = !h.arcs[2*h.nx].up;
			let firstarc = {
				let mut i = border-1;
				while i>0 && h.arcs[i].w<border {
					i -= 1
				}
				i
			};
			self.nx = h.nx - 1;
			for i in 1..firstarc {
				if h.arcs[i].w < firstarc {
					self.arcs[i] = h.arcs[i];
				} else if h.arcs[i].w < border {
					self.arcs[i] = h.arcs[i];
					self.arcs[i].w += 2*h.nx - border - 1;
					self.arcs[self.arcs[i].w] = Arc {
						w: i,
						up: self.arcs[i].up,
					};
				} else {
					self.arcs[i] = h.arcs[i];
					self.arcs[i].w -= (h.arcs[firstarc].w-border) + (border-firstarc);
					self.arcs[self.arcs[i].w] = Arc {
						w: i,
						up: self.arcs[i].up,
					};
				}
			}
			self.arcs[firstarc] = Arc {
				w: 2*h.nx - (border-firstarc) - 1,
				up: !side,
			};
			self.arcs[self.arcs[firstarc].w] = Arc {
				w: firstarc,
				up: !side,
			};
			for i in (firstarc+1)..border {
				if h.arcs[i].w < border && h.arcs[i].w > i {
					let a1 = i+2*h.nx-border-1;
					let a2 = h.arcs[i].w+2*h.nx-border-1;
					self.arcs[a1] = Arc {
						w: a2,
						up: h.arcs[i].up,
					};
					self.arcs[a2] = Arc {
						w: a1,
						up: self.arcs[a1].up,
					};
				}
			}
			for i in (border+1)..(2*h.nx) {
				if h.arcs[i].w > border && i < h.arcs[i].w {
					let a1 = if i < h.arcs[firstarc].w {
						2*h.nx-(h.arcs[firstarc].w-i)-(border-firstarc)-1
					} else {
						firstarc+(i-h.arcs[firstarc].w)
					};
					let j = h.arcs[i].w;
					let a2 = if j < h.arcs[firstarc].w {
						2*h.nx-(h.arcs[firstarc].w-j)-(border-firstarc)-1
					} else {
						firstarc+(j-h.arcs[firstarc].w)
					};
					self.arcs[a1] = Arc {
						w: a2,
						up: h.arcs[i].up,
					};
					self.arcs[a2] = Arc {
						w: a1,
						up: self.arcs[a1].up,
					};
				}
			}
			(ind[border] - 1, side)
		}
	}
	
	fn normalize_curve(&mut self, ind: &mut Numeration) {
		if self.nx <= 2 { return };
		let (index, side) = self.kill_arc(ind);
		debug_assert_ne!(index, 0);
		self.put_indices(ind);
		self.normalize_curve(ind);
		let h = self.clone();
		let i = ind.iter().position(|&x| x == index).unwrap();
		self.make_arc(&h, i, side);
		self.put_indices(ind);
	}
	
	fn overturn_curve(&mut self) {
		let h = self.clone();
		let n = 2*self.nx+1;
		for i in 1..n {
			let j = n-i;
			self.arcs[i] = Arc {
				w: n-h.arcs[j].w,
				up: h.arcs[j].up
			};
		}
		self.dir_r = !h.dir_r;
	}
	
	fn _gauss(&self) -> Vec<usize> {
		let mut v = vec!(0);
		let mut i = 1;
		for _ in 1..=self.nx {
			let j = self.arcs[i].w;
			v.push(i.min(j));
			i = j - 1;
			let j = self.arcs[i].w;
			v.push(i.min(j));
			i = j + 1;
		}
		v
	}
	
	fn _gauss_paired(&self) -> Vec<usize> {
		let mut g = self._gauss();
		let mut indices: Vec<Option<usize>> = vec![None; self.nx*2+1];
		for i in 1..=(2*self.nx) {
			let x = g[i];
			match indices[x] {
				None => indices[x] = Some(i),
				Some(j) => (g[i], g[j]) = (j, i),
			}
		}
		g
	}
	
	fn nonint(&self) -> bool {
		let g = self.arcs; // could use self._gauss_paired() as well
		for ifrom in 1..(2*self.nx) {
			let ito = g[ifrom].w;
			for j in (ifrom+1)..=(ito-1) {
				let jto = g[j].w;
				if jto < ifrom || jto > ito {
					return false
				}
			}
		}
		true
	}
/* In Pascal:
function nonint(c: curve): boolean;
var
  ifrom, ito, j, jto: longint;
  intersection: boolean;
begin
  intersection := false;
  for ifrom := 1 to 2*c.NumCross do
  begin
    ito := c.arcs[ifrom].where;
    if ito-1 < ifrom+1 then continue;
    for j := ifrom+1 to ito-1 do
    begin
      jto := c.arcs[j].where;
      if (jto < ifrom) or (jto > ito) then
      begin
        intersection := true;
        break;
      end;
    end;
    if intersection then break;
  end;
  nonint := not intersection;
end;
*/
	fn index(&self) -> i32 {
		let mut ind = 0;
		let mut i = 1;
		for _ in 1..=self.nx {
			let j = self.arcs[i].w;
			ind += if (j > i) ^ self.arcs[i].up {1} else {-1};
			i = j - 1;
			let j = self.arcs[i].w;
			ind += if (j > i) ^ self.arcs[i].up {1} else {-1};
			i = j + 1;
		}
		debug_assert_eq!(i, 2*self.nx+1);
		debug_assert_eq!(ind%2, 0);
		(ind/2 + 1) * if self.dir_r {1} else {-1}
	}
}

#[derive(Debug, PartialEq, Eq)]
struct DataCounts { c: [[usize; 4*MAX_X+1]; MAX_X+1]}

impl From<[[usize; 4*MAX_X+1]; MAX_X+1]> for DataCounts {
	fn from(c: [[usize; 4*MAX_X+1]; MAX_X+1]) -> Self {
		Self { c }
	}
}

impl DataCounts {
	fn check_data(&self, max: usize) {
		for i in 0..=max {
			for j in 1..=(4*max) {
				assert_eq!(self.c[i][j] % j, 0);
			}
		}
	}

	fn count_totals(&self, max: usize, divide: bool) -> Vec<usize> {
		self.check_data(max);
		(0..=max).map(|i| (1..=(4*max)).map(|j| self.c[i][j] / if divide {j} else {1}).sum()).collect()
	}
}

#[derive(Debug, PartialEq, Eq)]
struct Counts {
	pc: DataCounts,
	pn: DataCounts,
	nc: DataCounts,
	nn: DataCounts,
	by_index: [[[usize; MAX_X+2]; 4*MAX_X+1]; MAX_X+1],
}

impl Counts {
	fn new() -> Self {
		let mut s = Self {
			pc: [[0; 4*MAX_X+1]; MAX_X+1].into(),
			pn: [[0; 4*MAX_X+1]; MAX_X+1].into(),
			nc: [[0; 4*MAX_X+1]; MAX_X+1].into(),
			nn: [[0; 4*MAX_X+1]; MAX_X+1].into(),
			by_index: [[[0; MAX_X+2]; 4*MAX_X+1]; MAX_X+1],
		};
		s.pc.c[0][1] = 2;
		s.pn.c[0][1] = 1;
		s.nc.c[0][1] = 1;
		s.nn.c[0][1] = 1;
		s.by_index[0][1][1] = 1;
		s
	}
	
	fn print(&self, max: usize) {
		println!("Closed curves:");
		println!(" Plane is oriented, circle is oriented (OO):");
		println!("{:?}", self.pc.count_totals(max, true));
		println!(" Only plane is oriented:");
		println!("{:?}", self.pn.count_totals(max, true));
		println!(" Only circle is oriented:");
		println!("{:?}", self.nc.count_totals(max, true));
		println!(" Nothing is oriented:");
		println!("{:?}", self.nn.count_totals(max, true));
		println!(" OO by (nonnegative) index:");
		for i in 0..=max {
			let _ = (0..=(i+1)).map(|ind| (1..=(4*max)).map(move |j| assert_eq!(self.by_index[i][j][ind] % j, 0)));
			let c: Vec<usize> = (0..=(i+1)).map(|ind| (1..=(4*max)).map(|j| self.by_index[i][j][ind] / j).sum()).collect();
			println!("{c:?}");
		}
		println!("Long curves:");
		println!("{:?}", self.nn.count_totals(max, false));
	}
	
	fn inc(&mut self, current: usize, i: usize, f: bool, coef: usize, coefs: u8, index: usize) {
		self.by_index[current][i][index] += if index == 0 && f {2} else {1};
		if !f {
			self.pc.c[current][i] += 1;
			self.pn.c[current][i] += 1;
			self.nc.c[current][coef*i] += 1;
			self.nn.c[current][coef*i] += 1;
		} else if coefs == 1 {
            self.pc.c[current][i] += 2;
            self.pn.c[current][i] += 1;
            self.nc.c[current][i] += 1;
            self.nn.c[current][i] += 1;
		} else {
            self.pc.c[current][i] += 2;
            self.pn.c[current][i] += 1;
            self.nc.c[current][i] += [0, 2, 1][coef];
            self.nn.c[current][coef*i] += 1;
		}
	}
}

#[derive(Debug)]
struct Data {
	max: usize,
	curves: [Curve; MAX_X+1],
	indices: [Numeration; MAX_X+1],
	up_connect: [Free; MAX_X+1],
	down_connect: [Free; MAX_X+1],
	counts: Counts,
	counts_nonintersecting: Counts,
	current: usize,
}

impl Data {
	fn new(max: usize) -> Self {
		assert!(max > 0 && max <= MAX_X);
		let mut d = Self {
			max,
			curves: core::array::from_fn(|_| Curve::default()),
			indices: [[0; 2*MAX_X + 2]; MAX_X+1],
			up_connect: [[false; 2*MAX_X + 2]; MAX_X+1],
			down_connect: [[false; 2*MAX_X + 2]; MAX_X+1],
			counts: Counts::new(),
			counts_nonintersecting: Counts::new(),
			current: 0,
		};
		d.curves[0].put_indices(&mut d.indices[0]);
		d.curves[0].look_for_free_segments(&mut d.up_connect[1], true);
		d.curves[0].look_for_free_segments(&mut d.down_connect[1], false);
		d
	}
	
	fn curve_up(&mut self, segment: usize, side: bool) {
		self.current += 1;
		if self.current < self.max {
			self.curves[self.current+1].nx = 0;
		}
		let (c1, c2) = self.curves.split_at_mut(self.current);
		let cc = &mut c2[0];
		cc.make_arc(&c1[self.current-1], segment, side);
		cc.put_indices(&mut self.indices[self.current]);
		let mut h = cc.clone();
		let mut ind = self.indices[self.current];
		let mut symm = h.clone();
		symm.overturn_curve();
		let mut symind: Numeration = [0; 2*MAX_X + 2];
		symm.put_indices(&mut symind);
		symm.normalize_curve(&mut symind);
		let mut i = 0;
		let (mut coef, mut coefs) = (2, 2);
		while i == 0 || h != *cc {
			h.remake_curve(&ind);
			h.put_indices(&mut ind);
			h.normalize_curve(&mut ind);
			if coef == 2 && h == symm {
				coef = 1;
				if h.dir_r != symm.dir_r {
					coefs = 1;
				}
			}
			i += 1;
		}
		let flag = h.dir_r == cc.dir_r;
		let index = cc.index();
		if h.nonint() {
			self.counts_nonintersecting.inc(self.current, i, flag, coef, coefs, index.abs() as usize);
		}
		self.counts.inc(self.current, i, flag, coef, coefs, index.abs() as usize);
		cc.look_for_free_segments(&mut self.up_connect[self.current], true);
		cc.look_for_free_segments(&mut self.down_connect[self.current], false);
	}
	
	fn count(&mut self) {
		self.curve_up(1, true);
		loop {
			if self.current == self.max {
				self.current -= 1;
				continue
			}
			if self.curves[self.current+1].nx == 0 {
				self.curve_up(2*self.current+1, true);
				continue
			}
			let mut i = self.curves[self.current+1].arcs[2*self.current+2].w-1;
			if self.curves[self.current+1].arcs[2*self.current+2].up {
				while i > 0 && !(self.up_connect[self.current][i] && self.indices[self.current][i]%2 == 1) {
					i -= 1;
				}
				if i == 0 {
					self.curve_up(2*self.current+1, false);
				} else {
					self.curve_up(i, true);
				}
			} else {
				while i > 0 && !(self.down_connect[self.current][i] && self.indices[self.current][i]%2 == 1) {
					i -= 1;
				}
				if i == 0 {
					if self.current == 0 {
						break
					}
					self.current -= 1;
				} else {
					self.curve_up(i, false);
				}
			}
		}
	}
	
	fn print(&self) {
		println!("== With intersections in the Gauss diagram ==");
		self.counts.print(self.max);
		println!("== Without intersections in the Gauss diagram ==");
		self.counts_nonintersecting.print(self.max);
	}
}

fn main() {
	let max = std::env::args().skip(1).next().map_or(6, |a| a.parse::<usize>().expect("Command line argument must be an integer"));
	let mut d = Data::new(max);
	d.count();
    d.print();
}

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn test() {
		let mut d = Data::new(6);
		d.count();
		let Counts { pc, pn, nc, nn, .. } = d.counts;
		assert_eq!(&pc.count_totals(d.max, true), &[2, 3, 10, 39, 204, 1262, 8984], "A008980");
		assert_eq!(&pn.count_totals(d.max, true), &[1, 2, 5, 21, 102, 639, 4492], "A008981");
		assert_eq!(&nc.count_totals(d.max, true), &[1, 2, 5, 21, 102, 640, 4492], "A008982");
		assert_eq!(&nn.count_totals(d.max, true), &[1, 2, 5, 20, 82, 435, 2645], "A008983");
		assert_eq!(&nn.count_totals(d.max, false), &[1, 2, 8, 42, 260, 1796, 13396], "A054993");
		assert_eq!(&d.counts_nonintersecting.nn.count_totals(4, true), &[1, 2, 5, 18, 70], "A118814");
	}
}
