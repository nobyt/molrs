//! 2D 点/ベクトル (D1)。
//!
//! レイアウトは結合長 = 1.0 の無次元単位で行う。角度はラジアン、
//! x 軸正方向 = 0、反時計回り正 (数学的慣例。SVG 出力時に y を反転する)。

use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    pub const ZERO: Point2 = Point2 { x: 0.0, y: 0.0 };

    pub fn new(x: f64, y: f64) -> Point2 {
        Point2 { x, y }
    }

    /// 角度 theta (rad) 方向の単位ベクトル。
    pub fn from_angle(theta: f64) -> Point2 {
        Point2::new(theta.cos(), theta.sin())
    }

    pub fn dot(self, o: Point2) -> f64 {
        self.x * o.x + self.y * o.y
    }

    /// 2D 外積 (z 成分)。正 = o が self の反時計回り側。
    pub fn cross(self, o: Point2) -> f64 {
        self.x * o.y - self.y * o.x
    }

    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    pub fn distance(self, o: Point2) -> f64 {
        (self - o).norm()
    }

    pub fn normalized(self) -> Option<Point2> {
        let n = self.norm();
        if n < 1e-12 {
            None
        } else {
            Some(self / n)
        }
    }

    /// 反時計回りに 90° 回した垂直ベクトル。
    pub fn perp(self) -> Point2 {
        Point2::new(-self.y, self.x)
    }

    /// x 軸からの角度 (rad, (-π, π])。
    pub fn angle(self) -> f64 {
        self.y.atan2(self.x)
    }

    /// 原点まわりに theta (rad) 回転。
    pub fn rotated(self, theta: f64) -> Point2 {
        let (s, c) = theta.sin_cos();
        Point2::new(c * self.x - s * self.y, s * self.x + c * self.y)
    }

    /// 中心 c まわりに theta (rad) 回転。
    pub fn rotated_about(self, c: Point2, theta: f64) -> Point2 {
        (self - c).rotated(theta) + c
    }
}

impl Add for Point2 {
    type Output = Point2;
    fn add(self, o: Point2) -> Point2 {
        Point2::new(self.x + o.x, self.y + o.y)
    }
}

impl Sub for Point2 {
    type Output = Point2;
    fn sub(self, o: Point2) -> Point2 {
        Point2::new(self.x - o.x, self.y - o.y)
    }
}

impl Neg for Point2 {
    type Output = Point2;
    fn neg(self) -> Point2 {
        Point2::new(-self.x, -self.y)
    }
}

impl Mul<f64> for Point2 {
    type Output = Point2;
    fn mul(self, k: f64) -> Point2 {
        Point2::new(self.x * k, self.y * k)
    }
}

impl Div<f64> for Point2 {
    type Output = Point2;
    fn div(self, k: f64) -> Point2 {
        Point2::new(self.x / k, self.y / k)
    }
}

/// 角度を最も近い 30° の倍数 (rad) にスナップする。
/// IUPAC 2008 は「水平から ±30° の結合を最大化する」配向を推奨しており、
/// 鎖結合の量子化に使う。
pub fn snap_angle_30(theta: f64) -> f64 {
    let step = std::f64::consts::PI / 6.0;
    (theta / step).round() * step
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_6, PI};

    const EPS: f64 = 1e-12;

    #[test]
    fn arithmetic() {
        let a = Point2::new(1.0, 2.0);
        let b = Point2::new(3.0, -1.0);
        assert_eq!(a + b, Point2::new(4.0, 1.0));
        assert_eq!(a - b, Point2::new(-2.0, 3.0));
        assert_eq!(-a, Point2::new(-1.0, -2.0));
        assert_eq!(a * 2.0, Point2::new(2.0, 4.0));
        assert_eq!(a / 2.0, Point2::new(0.5, 1.0));
        assert!((a.dot(b) - 1.0).abs() < EPS);
        assert!((a.cross(b) + 7.0).abs() < EPS);
    }

    #[test]
    fn rotation() {
        let e = Point2::new(1.0, 0.0);
        let r = e.rotated(FRAC_PI_2);
        assert!((r.x).abs() < EPS && (r.y - 1.0).abs() < EPS);
        // 120° 回転 3 回で一周
        let mut p = Point2::new(1.0, 0.5);
        for _ in 0..3 {
            p = p.rotated(2.0 * FRAC_PI_3);
        }
        assert!(p.distance(Point2::new(1.0, 0.5)) < 1e-9);
        // 中心指定回転
        let q = Point2::new(2.0, 0.0).rotated_about(Point2::new(1.0, 0.0), PI);
        assert!(q.distance(Point2::ZERO) < EPS);
    }

    #[test]
    fn angles_and_snap() {
        assert!((Point2::new(0.0, 1.0).angle() - FRAC_PI_2).abs() < EPS);
        assert!((Point2::from_angle(FRAC_PI_6).angle() - FRAC_PI_6).abs() < EPS);
        // 29° → 30°, 44° → 30°? (44° は 30° と 60° の中間より 60° 側)
        let deg = |d: f64| d * PI / 180.0;
        assert!((snap_angle_30(deg(29.0)) - deg(30.0)).abs() < EPS);
        assert!((snap_angle_30(deg(46.0)) - deg(60.0)).abs() < EPS);
        assert!((snap_angle_30(deg(-14.0)) - 0.0).abs() < EPS);
        assert!((snap_angle_30(deg(-16.0)) - deg(-30.0)).abs() < EPS);
    }

    #[test]
    fn perp_and_norm() {
        let v = Point2::new(3.0, 4.0);
        assert!((v.norm() - 5.0).abs() < EPS);
        assert!((v.perp().dot(v)).abs() < EPS);
        assert!(v.perp().cross(v) < 0.0);
        assert!(v.normalized().unwrap().norm() - 1.0 < EPS);
        assert!(Point2::ZERO.normalized().is_none());
    }
}
